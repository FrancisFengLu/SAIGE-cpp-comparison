#define ARMA_USE_SUPERLU 1
#include <RcppArmadillo.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <mutex>
#include <cmath>
#include <ctime>// include this header for calculating execution time
#include <cassert>
#include <random>  // for std::random_device, std::mt19937, std::bernoulli_distribution
#include <boost/date_time.hpp> // for gettimeofday and timeval
#include "getMem.hpp"
#include "UTIL.hpp"  // Substituted from src/UTIL.hpp for utility functions
#include "SAIGE_step1_fast.hpp"  // Included from src/Main.hpp for function declarations
#include "bed_reader.hpp"      // PR-5: parallel pread-backed BED reader
#include "marker_decoder.hpp"  // PR-5: decode + QC + repack
#include "parallel_decode.hpp" // PR-5: parallel_decode_bed()
#include "packed_store.hpp"    // option 3: PackedFlat primary storage
#include "gpu_matvec.hpp"      // G3: optional cuBLAS K·u acceleration
#include "tools/avx2_kernel/avx2_kernel.hpp"  // Phase-1: fused 2-bit decode kernels
#include <RcppParallel.h>
#include <RcppParallel/TBB.h>
#include <cstdlib>

// Optional R-comparison bypass reads / debug dumps. Opt-in via environment
// variables (SAIGE_BYPASS_DIR / SAIGE_DEBUG_DIR); returns "" when unset so the
// gated readers/writers below become no-ops. Default: off (nothing written/read).
static inline std::string saige_env_path(const char* var, const std::string& name) {
  const char* e = std::getenv(var);
  return (e && *e) ? (std::string(e) + "/" + name) : std::string();
}
using namespace Rcpp;
using namespace std;
using namespace RcppParallel;


// R CONNECTION: Global variable used in R functions for quality control thresholds
float minMAFtoConstructGRM = 0;
// R CONNECTION: This C++ class stores genotype data and is initialized via setgeno() 
// which is called from R functions like SAIGE_fitNULLGLMM() and SAIGE_fitNULLGLMM_fast()
// The class methods are accessed through various Rcpp::export functions
//This is a class with attritbutes about the genotype informaiton 
class genoClass{
private:
        //COPY from RVTEST:
        // we reverse the two bits as defined in PLINK format
        const static unsigned char HOM_REF = 0x0;  // 0b00 ;
        const static unsigned char HET = 0x2;      // 0b10 ;
        const static unsigned char HOM_ALT = 0x3;  // 0b11 ;
        const static unsigned char MISSING = 0x1;  // 0b01 ;


public:
        //to chunk the geno vector to avoid large continuous memory usage 
	// R CONNECTION: Memory management parameters set by R based on memoryChunk argument
	int numMarkersofEachArray;
        int numofGenoArray;
        int numMarkersofLastArray;
        std::vector< std::vector<unsigned char>* > genoVecofPointers;  // legacy storage, only populated by the serial BED loader (SAIGE_SERIAL_BED=1)
        ///////////
        std::vector< std::vector<unsigned char>* > genoVecofPointers_forVarRatio;  // ditto, VR pool

        // Option-3 primary storage for pass-QC packed bytes. Populated in the
        // non-VR parallel path (setGenoObj) via std::move from parallel_decode_bed().
        // When use_packed_flat_ is true, readers use this instead of
        // genoVecofPointers[idx]->at(j), which avoids the 2× peak RSS seen in PR-5.
        saige::PackedFlat packed_flat_;
        bool              use_packed_flat_       = false;
        bool              packed_flat_released_  = false;  // G4: host copy freed after GPU upload

        // Same storage trick for the (tiny, ≤1000-marker) variance-ratio pool.
        // Populated by the parallel BED path; replaces the per-marker
        // genoVecofPointers_forVarRatio heap vectors the serial path builds.
        saige::PackedFlat packed_flat_vr_;
        bool              use_packed_flat_vr_    = false;

        // Return the j-th packed byte of the snp-th pass-QC marker, routing
        // between the option-3 flat buffer and the legacy vector-of-pointers.
        inline unsigned char packed_byte(std::size_t snp_idx, std::size_t byte_idx) const {
          if (packed_flat_released_) {
            // G4: we freed the host copy when GPU took ownership of the data.
            // Any caller reaching this point wants a host-side byte, which we
            // no longer have. Surfacing a clear error beats silently reading
            // freed memory.
            throw std::runtime_error(
                "genoClass::packed_byte called after packed_flat_ was released "
                "to the GPU backend (--gpu mode). No CPU fallback possible; "
                "re-run without --gpu or disable the release hook.");
          }
          if (use_packed_flat_) {
            return packed_flat_.raw()[snp_idx * packed_flat_.nbyte() + byte_idx];
          }
          return genoVecofPointers[snp_idx]->at(byte_idx);
        }
        // Same routing for the variance-ratio pool. Never released to the GPU
        // (VR estimation is host-side), so no released_ guard is needed.
        inline unsigned char packed_byte_vr(std::size_t snp_idx, std::size_t byte_idx) const {
          if (use_packed_flat_vr_) {
            return packed_flat_vr_.raw()[snp_idx * packed_flat_vr_.nbyte() + byte_idx];
          }
          return genoVecofPointers_forVarRatio[snp_idx]->at(byte_idx);
        }
        inline std::size_t packed_size(std::size_t snp_idx) const {
          if (packed_flat_released_) return packed_flat_.nbyte();  // size is still known, bytes aren't
          if (use_packed_flat_) return packed_flat_.nbyte();
          return genoVecofPointers[snp_idx]->size();
        }
        // Number of addressable marker rows in whichever store is live.
        inline std::size_t packed_n_markers() const {
          if (use_packed_flat_) return packed_flat_.n_stored();
          return genoVecofPointers.size();
        }

        // Phase-1 AVX2 kernel support: contiguous pointer to the packed bytes
        // of one pass-QC marker. Valid in both storage modes because
        // numMarkersofEachArray == 1 (each legacy vector holds exactly one
        // marker's ⌈N/4⌉ bytes). Callers must check packed_rows_contiguous().
        inline const unsigned char* packed_row_ptr(std::size_t snp_idx) const {
          if (use_packed_flat_) {
            return packed_flat_.raw() + snp_idx * packed_flat_.nbyte();
          }
          return genoVecofPointers[snp_idx]->data();
        }
        inline bool packed_rows_contiguous() const {
          if (packed_flat_released_) return false;
          if (use_packed_flat_) return true;
          return numMarkersofEachArray == 1 && !genoVecofPointers.empty();
        }

        // G4: release the host-side packed bytes once the GPU backend has
        // copied them to device memory. Frees ~6 GB on UKB 100K. After this
        // call, packed_byte() throws — the caller is committed to GPU for
        // the rest of the run.
        void release_packed_flat_host() {
          if (packed_flat_released_) return;
          // Move-assign an empty PackedFlat over it: the old buffer is freed
          // by the assignment. The shell keeps the same nbyte() so the
          // packed_size() call sites that only want the geometry still work.
          {
            saige::PackedFlat tmp_shell;
            // Preserve the reported nbyte() for packed_size() queries even
            // after the data is gone. PackedFlat doesn't expose a setter, so
            // we give tmp_shell the same geometry without allocating.
            tmp_shell.init(/*capacity_markers=*/0, packed_flat_.nbyte());
            packed_flat_ = std::move(tmp_shell);  // old packed_flat_ destructs, frees memory
          }
          packed_flat_released_ = true;
          std::cout << "[genoClass] packed_flat_ host buffer released "
                    << "(GPU mode — ~"
                    << (static_cast<std::size_t>(packed_flat_.nbyte()) *
                        static_cast<std::size_t>(use_packed_flat_ ? 1 : 0)) / (1024ULL*1024ULL)
                    << " MB freed)\n";
        }
	//arma::fvec g_cateVarRatioMinMACVecExclude;
	//arma::fvec g_cateVarRatioMaxMACVecInclude;
	float g_minMACVarRatio;
	float g_maxMACVarRatio;
	bool isVarRatio = false;
	int numberofMarkers_varRatio = 0;
	int numberofMarkers_varRatio_common = 0;
	arma::ivec g_randMarkerIndforVR;
	std::vector<float>      invstdvVec0_forVarRatio;
        arma::fvec      invstdvVec_forVarRatio;
	 std::vector<float>      alleleFreqVec0_forVarRatio;
        arma::fvec      alleleFreqVec_forVarRatio;
	std::vector<int>      MACVec0_forVarRatio;
	std::vector<int>      markerIndexVec0_forVarRatio;
	arma::ivec MACVec_forVarRatio;
	arma::ivec markerIndexVec_forVarRatio;


	//vector<unsigned char> genoVec; 	 
  	size_t M;
  	size_t N;
	size_t Nnomissing;
	// R CONNECTION: Inverse standard deviation vectors used for genotype standardization, accessed by R via getAlleleFreqVec()
	std::vector<float>	invstdvVec0;
	arma::fvec	invstdvVec;
	vector<int>	ptrsubSampleInGeno;
	std::vector<bool> indicatorGenoSamplesWithPheno_in;	
	

  	// R CONNECTION: Allele frequency vectors returned to R via getAlleleFreqVec() for QC and analysis
  	std::vector<float> 	alleleFreqVec0;
  	arma::fvec 	alleleFreqVec;
  	arma::ivec	m_OneSNP_Geno;
  	arma::fvec	m_OneSNP_StdGeno;
  	arma::fvec	m_DiagStd;
	arma::fvec	m_DiagStd_LOCO;
  	arma::fmat	mtx_DiagStd_LOCO;


	std::vector<int>	MACVec0; //for variance ratio based on different MAC categories
	arma::ivec	MACVec;
	std::vector<int>	origPlinkIdx0;  // mapping: main array index -> original plink marker index
	arma::ivec	subMarkerIndex; //for sparse GRM
	arma::fmat      stdGenoMultiMarkersMat;	
	std::vector<float> stdGenoforSamples; //for sparse GRM
	std::vector<float>     kinValueVecFinal;
        float relatednessCutoff;
	float maxMissingRate;

	std::vector< std::pair<int, int> > indiceVec;
	std::vector<float> kinValueVecSparse;  // kinship values for sparse GRM pairs
	arma::ivec xout;
        arma::ivec yout;
	//int Mmafge1perc;
	bool setKinDiagtoOne;
	int numberofMarkerswithMAFge_minMAFtoConstructGRM = 0;
//	arma::SpMat<float> sparseGRMinC(2,2);
	std::vector<bool> MarkerswithMAFge_minMAFtoConstructGRM_indVec;	


        //std::vector<float> stdGenoVec;
	//for LOCO
	//bool LOCO = false;
	//vector<int> chromosomeStartIndex;
	//vector<int> chromosomeEndIndex;
	//vector<int> chromosomeVec;
        size_t Msub;
        int startIndex;
        int endIndex;
	int chromIndex;

        
        arma::ivec startIndexVec;
        arma::ivec endIndexVec;
        arma::ivec startIndexVec_forvr;
        arma::ivec endIndexVec_forvr;


        int Msub_MAFge_minMAFtoConstructGRM;

	int Msub_MAFge_minMAFtoConstructGRM_singleChr;
	arma::ivec Msub_MAFge_minMAFtoConstructGRM_byChr;
	//end for LOCO

	unsigned char m_genotype_buffer[4];
	int geno_idx;
	int m_size_of_esi;
	unsigned char m_bits_val[8];

	
	//look-up table for std geno
	//float stdGenoLookUpArr[3] = {0};
	void setStdGenoLookUpArr(float mafVal, float invsdVal, arma::fvec & stdGenoLookUpArr){
	//	arma::fvec stdGenoLookUpArr(3);
		float mafVal2 = 2*mafVal;
		stdGenoLookUpArr(0) = (0-mafVal2)*invsdVal;
		stdGenoLookUpArr(1) = (1-mafVal2)*invsdVal;
		stdGenoLookUpArr(2) = (2-mafVal2)*invsdVal;
	//	return(stdGenoLookUpArr)
	}


        //look-up table in a 2D array for sparseKin 
        float sKinLookUpArr[3][3] = {{0}};
	//(g - 2*freq)* invStd;;
        void setSparseKinLookUpArr(float mafVal, float invsdVal){
		float mafVal2 = 2*mafVal;
		float a0 = (0-mafVal2)*invsdVal;
		float a1 = (1-mafVal2)*invsdVal;
		float a2 = (2-mafVal2)*invsdVal;
		
		sKinLookUpArr[0][0] = a0*a0;
		sKinLookUpArr[0][1] = a0*a1;
		sKinLookUpArr[0][2] = a0*a2;
		sKinLookUpArr[1][0] = sKinLookUpArr[0][1];
		sKinLookUpArr[1][1] = a1*a1;
		sKinLookUpArr[1][2] = a1*a2;
		sKinLookUpArr[2][0] = sKinLookUpArr[0][2];
		sKinLookUpArr[2][1] = sKinLookUpArr[1][2];
		sKinLookUpArr[2][2] = a2*a2;

	}




        void setBit(unsigned char & ch, int ii, int aVal, int bVal){

                if (bVal == 1 && aVal == 1){
			ch ^= char(1 << ((ii*2) + 1)); //set a to be 1

                }else if(bVal == 0){
			ch ^= char(1 << (ii*2)); //change b to 0

                        if(aVal == 1){
				ch ^= char(1 << ((ii*2) + 1)); //change a to 1
                        }
                }
        }



	//COPY from RVTEST:
	void setGenotype(unsigned char* c, const int pos, const int geno) {
    		(*c) |= (geno << (pos << 1));
  	}

	void getGenotype(unsigned char* c, const int pos, int& geno) {
    		geno = ((*c) >> (pos << 1)) & 0x3;  // 0b11 = 0x3
  	}



	void Init_OneSNP_Geno(){
		m_size_of_esi = (Nnomissing+3)/4;
		int k = 8;
		while (k > 0){
			-- k;
			m_bits_val[k] = 1 << k;
		}
	}
	

        arma::ivec * Get_OneSNP_Geno(size_t SNPIdx){
                m_OneSNP_Geno.zeros(Nnomissing);

		//avoid large continuous memory usage
		int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
		////////////////

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind= 0;
                unsigned char geno1;
                int bufferGeno;
                for(size_t i=Start_idx; i< Start_idx+m_size_of_esi - 1; i++){
                        //geno1 = genoVec[i];
			geno1 = packed_byte(indexOfVectorPointer, i); //option-3: flat or legacy
                        for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
				bufferGeno = 2-(a+b);
				m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                //if(ind >= Nnomissing){

                                ////printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                //        return & m_OneSNP_Geno;
                                //}
                        }
                }

		size_t i = Start_idx+m_size_of_esi - 1;
		geno1 = packed_byte(indexOfVectorPointer, i);
		for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                bufferGeno = 2-(a+b);
                                m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                if(ind >= Nnomissing){

                                ////printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                        return & m_OneSNP_Geno;
                                }
                }

                return & m_OneSNP_Geno;
       }
   
        arma::ivec * Get_OneSNP_Geno_forVarRatio(size_t SNPIdx){
                m_OneSNP_Geno.zeros(Nnomissing);

		//avoid large continuous memory usage
		int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
		////////////////

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind= 0;
                unsigned char geno1;
                int bufferGeno;
                for(size_t i=Start_idx; i< Start_idx+m_size_of_esi-1; i++){
                        //geno1 = genoVec[i];
			geno1 = packed_byte_vr(indexOfVectorPointer, i);
                        for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
				bufferGeno = 2-(a+b);
				m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                //if(ind >= Nnomissing){

                                //printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                //        return & m_OneSNP_Geno;
                                //}
                        }
                }

		size_t i = Start_idx+m_size_of_esi-1;
		geno1 = packed_byte_vr(indexOfVectorPointer, i);
                for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                bufferGeno = 2-(a+b);
                                m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                if(ind >= Nnomissing){

                                //printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                        return & m_OneSNP_Geno;
                                }
                  }

                return & m_OneSNP_Geno;
       }


	// Precomputed BED lookup table (built in setGenoObj)
	int (*bed_lookup_ptr)[4] = nullptr;
	void set_bed_lookup(int (*lut)[4]) { bed_lookup_ptr = lut; }

	void Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld, float & altFreq, float & missingRate, int & mac,  int & alleleCount, bool & passQC, size_t SNPIdx_new, bool & passVarRatio , size_t SNPIdx_vr){

		// Decode BED bytes directly, only track phenotyped samples
		// Avoids allocating N-length temp vector (N can be 408K)
		m_OneSNP_Geno.zeros(Nnomissing);
		const int nbytesFull = (N+3)/4;
		alleleCount = 0;
		int numMissing = 0;
		size_t ind = 0;  // FAM sample index

		// Temporary storage for all N genotypes (needed for subsetting later)
		// Use thread-local static to avoid per-call allocation
		static thread_local std::vector<int8_t> genoAll;
		if (genoAll.size() < N) genoAll.resize(N);

		// Decode all bytes using lookup table
		for (int i = 0; i < nbytesFull; i++) {
			const unsigned char byte = genoVecOneMarkerOld[i];
			const int* lut = bed_lookup_ptr[byte];
			const size_t base = (size_t)i * 4;
			for (int j = 0; j < 4 && (base + j) < N; j++) {
				const size_t fam_idx = base + j;
				const int geno_val = lut[j];
				genoAll[fam_idx] = geno_val;
				if (indicatorGenoSamplesWithPheno_in[fam_idx]) {
					if (geno_val == 3) {
						numMissing++;
					} else {
						alleleCount += geno_val;
					}
				}
			}
		}


	      altFreq = alleleCount/float((Nnomissing-numMissing) * 2);
	      //sum = 0;
	      //std::cout << "missingRate " << missingRate << std::endl;
	      //std::cout << "maxMissingRate " << maxMissingRate << std::endl;
	      missingRate = numMissing/float(Nnomissing);	      

              //int indxInOut = 0;
	      //if(minMAFtoConstructGRM > 0){
              //if(altFreq >= minMAFtoConstructGRM && altFreq <= (1-minMAFtoConstructGRM) && missingRate <= maxMissingRate){
	      	int fillinMissingGeno = int(round(2*altFreq)); 
	
		if(numMissing > 0){
		       //for(int indx=0; indx < Nnomissing; indx++){
                                                //cout << "HERE5\n";
			//	u = indx & 3;
			//	bufferGeno = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
			alleleCount = alleleCount + fillinMissingGeno*numMissing;
			
			//		if(bufferGeno == 3){
			//			bufferGeno = fillinMissingGeno;
			//			alleleCount = alleleCount + bufferGeno;
			//		}
			//	}	
  			/*
				//setGenotype(&geno2, u, bufferGeno);
				if(bufferGeno == 0){
                                        setGenotype(&geno2, u, HOM_ALT);
                                }else if(bufferGeno == 1){
                                        setGenotype(&geno2, u, HET);
                                }else if(bufferGeno == 2){
                                        setGenotype(&geno2, u, HOM_REF);
                                }
				//else{
                                //        setGenotype(&geno1, u, MISSING);
                                        //m_OneSNP_Geno[j] = 0;  //12-18-2017
                                //}	

				if(u == 3 || indx == (Nnomissing-1)){
                                        genoVecofPointers[SNPIdx/numMarkersofEachArray]->push_back(geno2); //avoid large continuous memory usage
                                        geno2 = 0;
               			}
		*/	
		}
			//passQC = true;	
	     //}

	     altFreq = alleleCount/float(Nnomissing * 2);

	     unsigned char geno2 = 0;  // Must initialize to 0 since setGenotype uses OR
	     passQC = false;
	     passVarRatio = false;
	     float maf = std::min(altFreq, 1-altFreq);
	     mac = std::min(alleleCount, int(Nnomissing) * 2 - alleleCount);


		if(maf >= minMAFtoConstructGRM && missingRate <= maxMissingRate){
			passQC = true;
		}
		if(isVarRatio){
			if(g_maxMACVarRatio != -1){ //if estimating categorical variance ratios
			   if(mac >= g_minMACVarRatio && mac < g_maxMACVarRatio){
				passVarRatio = true;
				genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
				genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
			   }else if(mac >= g_maxMACVarRatio){
				   //randomly select 200 markers for estimating the variance ratio for the last MAC category	
				   //if(numberofMarkers_varRatio_common < 200){
				   	//if(static_cast<int>(SNPIdx) == 123){
					//	std::cout << "123" << std::endl;
					//	bool isIng_randMarkerIndforVR = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
					//	std::cout << "isIng_randMarkerIndforVR " << isIng_randMarkerIndforVR << std::endl;
					//}
				   	passVarRatio = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
					if(passVarRatio){
						//std::cout << "SNPIdx " << SNPIdx << std::endl;
						genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
						genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
				  		numberofMarkers_varRatio_common = numberofMarkers_varRatio_common + 1;
						//passVarRatio = false;
					}
				  //} 
			//	passVarRatio = true;	
			}
			}else{
				if(mac >= g_minMACVarRatio){
                                   //randomly select 200 markers for estimating the variance ratio for the last MAC category
                                   //if(numberofMarkers_varRatio_common < 200){
                                        //if(static_cast<int>(SNPIdx) == 123){
                                        //      std::cout << "123" << std::endl;
                                        //      bool isIng_randMarkerIndforVR = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
                                        //      std::cout << "isIng_randMarkerIndforVR " << isIng_randMarkerIndforVR << std::endl;
                                        //}
                                        passVarRatio = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
                                        if(passVarRatio){
                                                //std::cout << "SNPIdx " << SNPIdx << std::endl;
                                                genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
                                                genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                                                numberofMarkers_varRatio_common = numberofMarkers_varRatio_common + 1;
                                                //passVarRatio = false;
                                        }
                                  //}
                        //      passVarRatio = true;
                        }	
			

			}
			//avoid the overlap between markers for GRM and markers for variance ratio estimation	   
			if(passVarRatio){
				passQC = false;
			}
		      //}
		}

		if(passQC | passVarRatio){
			// Dump first marker's raw genotypes before re-encoding
			static bool dumped_first_marker = false;
			if (!dumped_first_marker && passQC) {
				dumped_first_marker = true;
				std::string dump_path = "/data/ukb_cpp_test/output/cpp_first_marker_raw.csv";
				std::ofstream df(dump_path);
				if (df.is_open()) {
					df << "pheno_idx,fam_idx,raw_geno_from_bed,fillinMissing\n";
					for (int d = 0; d < std::min((int)Nnomissing, 200); d++) {
						int fam_i = ptrsubSampleInGeno[d] - 1;
						int g = genoAll[fam_i];
						df << d << "," << fam_i << "," << g << "," << fillinMissingGeno << "\n";
					}
					df.close();
					std::cout << "[DUMP] First QC marker (SNPIdx=" << SNPIdx << " SNPIdx_new=" << SNPIdx_new
					          << ") raw genotypes: " << dump_path << std::endl;
					std::cout << "[DUMP] freq=" << altFreq << " N=" << Nnomissing
					          << " ptrsubSample[0:3]=" << ptrsubSampleInGeno[0] << "," << ptrsubSampleInGeno[1]
					          << "," << ptrsubSampleInGeno[2] << std::endl;
				}
			}
			int u;
			int bufferGeno;
			for(int indx=0; indx < Nnomissing; indx++){
                              u = indx & 3;
                              bufferGeno = genoAll[ptrsubSampleInGeno[indx] - 1];
			      if(bufferGeno == 3){
				bufferGeno = fillinMissingGeno;
			      }
			      if(bufferGeno == 0){
                                        setGenotype(&geno2, u, HOM_ALT);
                                }else if(bufferGeno == 1){
                                        setGenotype(&geno2, u, HET);
                                }else if(bufferGeno == 2){
                                        setGenotype(&geno2, u, HOM_REF);
                              }
			      if(u == 3 || indx == (Nnomissing-1)){
				       if(passVarRatio){
						genoVecofPointers_forVarRatio[SNPIdx_vr/numMarkersofEachArray]->push_back(geno2); //avoid large continuous memory usage
					}
				        if(passQC){
						genoVecofPointers[SNPIdx_new/numMarkersofEachArray]->push_back(geno2);
					}	
                                        geno2 = 0;
                              }
			}
		}
	    // altFreq = alleleCount/float(Nnomissing * 2);

   }
		

   int Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out ){
		//avoid large continuous memory usage
                int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
                ////////////////
		//std::cout << "indexOfVectorPointer " << indexOfVectorPointer << std::endl;
 		out->zeros(Nnomissing);
		//std::cout << "m_size_of_esi " << m_size_of_esi << std::endl;
		//std::cout << "SNPIdxinVec " << SNPIdxinVec << std::endl;
		//std::cout << "genoVecofPointers[indexOfVectorPointer]->size() " << genoVecofPointers[indexOfVectorPointer]->size() << std::endl;



 		size_t Start_idx = m_size_of_esi * SNPIdxinVec;

		//std::cout << "Start_idx " << Start_idx << std::endl;
		size_t ind= 0;
		unsigned char geno1;
		
		float freq = alleleFreqVec[SNPIdx];
		//cout << "Get_OneSNP_StdGeno here" << endl; 
		float invStd = invstdvVec[SNPIdx];

		arma::fvec stdGenoLookUpArr(3);
		setStdGenoLookUpArr(freq, invStd, stdGenoLookUpArr);
//		std::cout << "freq " << freq << endl;
//		std::cout << "invStd " << invStd << endl;

		//setStdGenoLookUpArr(freq, invStd);
		//std::cout << "stdGenoLookUpArr[0]: " << stdGenoLookUpArr[0] << std::endl;
		//std::cout << "stdGenoLookUpArr[1]: " << stdGenoLookUpArr[1] << std::endl;
		//std::cout << "stdGenoLookUpArr[2]: " << stdGenoLookUpArr[2] << std::endl;
//		cout << "Get_OneSNP_StdGeno here2"  << endl;
		for(size_t i=Start_idx; i< Start_idx+m_size_of_esi-1; i++){
//			geno1 = genoVec[i];
			geno1 = packed_byte(indexOfVectorPointer, i);

			for(int j=0; j<4; j++){
    			int b = geno1 & 1 ;
    			geno1 = geno1 >> 1;
    			int a = geno1 & 1 ;
    			//(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
    			(*out)[ind] = stdGenoLookUpArr(2-(a+b));
//			std::cout << "a " << a << endl;
//			std::cout << "b " << b << endl;
//			std::cout << "(*out)[ind] " << (*out)[ind] << endl;
			ind++;
    			geno1 = geno1 >> 1;
    			
//    			if(ind >= Nnomissing){
//				cout << "Get_OneSNP_StdGeno " << SNPIdx << endl; 
//				cout << "Nnomissing " << Nnomissing << endl; 
//				stdGenoLookUpArr.clear();
//    				return 1;
//    			}
	    		}
		}


		size_t i = Start_idx+m_size_of_esi-1;
                geno1 = packed_byte(indexOfVectorPointer, i);

                for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        (*out)[ind] = stdGenoLookUpArr(2-(a+b));
                        ind++;
                        geno1 = geno1 >> 1;

                        if(ind >= Nnomissing){
                                stdGenoLookUpArr.clear();
                                return 1;
                        }
                }
		//cout << "SNPIdx " << SNPIdx << endl; 

		stdGenoLookUpArr.clear();
		return 1;
				
 	}


	arma::fvec * Get_Diagof_StdGeno(){
	
		arma::fvec * temp = &m_OneSNP_StdGeno;
		// Not yet calculated
		//cout << "size(m_DiagStd)[0] " << size(m_DiagStd)[0] << endl;
		if(size(m_DiagStd)[0] != Nnomissing){
			m_DiagStd.zeros(Nnomissing);

			// DEBUG: Output stdGeno for first 3 markers to file
			std::ofstream stdgeno_file(saige_env_path("SAIGE_DEBUG_DIR", "cpp_stdgeno.txt"));
			stdgeno_file << "# StdGeno values for first 3 markers, all samples\n";
			stdgeno_file << "# Columns: Sample, Marker0, Marker1, Marker2\n";

			arma::fmat first3_stdgeno(Nnomissing, 3);

			// DEBUG: Track sample 0's stdGeno^2 accumulation
			float sample0_cumsum = 0;
			std::ofstream sample0_file(saige_env_path("SAIGE_DEBUG_DIR", "cpp_sample0_cumsum.txt"));
			sample0_file << "# Sample 0 stdGeno^2 cumulative sum by marker\n";
			sample0_file << "# Columns: Marker, stdGeno[0], stdGeno^2[0], cumsum\n";

			for(size_t i=0; i< numberofMarkerswithMAFge_minMAFtoConstructGRM; i++){
				//if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){


				Get_OneSNP_StdGeno(i, temp);

				// Save first 3 markers' stdGeno
				if(i < 3) {
					first3_stdgeno.col(i) = *temp;
				}

				/*if(i == 0){
					cout << "setgeno mark7 " << i <<  endl;
					for(int j=0; j<10; ++j)
					{
                				cout << (*temp)[j] << ' ';
                			}
                			cout << endl;
				}
				*/
				m_DiagStd = m_DiagStd + (*temp) % (*temp);

				// DEBUG: Track sample 0's cumsum
				float val0 = (*temp)[0];
				sample0_cumsum += val0 * val0;
				if(i < 100 || i % 1000 == 0) {
					sample0_file << i << "\t" << val0 << "\t" << val0*val0 << "\t" << sample0_cumsum << "\n";
				}

				//}
		//		std::cout << "i " << i << std::endl;
		//		std::cout << "numberofMarkerswithMAFge_minMAFtoConstructGRM " << numberofMarkerswithMAFge_minMAFtoConstructGRM << std::endl;
			}

			sample0_file << "FINAL\t-\t-\t" << sample0_cumsum << "\n";
			sample0_file.close();
			std::cout << "DEBUG: Sample 0 final cumsum = " << sample0_cumsum << std::endl;

			// Write stdGeno to file
			for(size_t j=0; j < Nnomissing; j++) {
				stdgeno_file << j << "\t" << first3_stdgeno(j,0) << "\t"
				             << first3_stdgeno(j,1) << "\t" << first3_stdgeno(j,2) << "\n";
			}
			stdgeno_file.close();
			std::cout << "DEBUG: Wrote stdGeno (SAIGE_DEBUG_DIR, if set)" << std::endl;

		}
/*
		std::cout << "test\n";
		for(int i=0; i<10; ++i)
        	{
        	  cout << m_DiagStd[i] << ' ';
        	}
		cout << endl;
*/
		// DEBUG: Override sample 0's diagonal with R's value to test tau matching
		// R's normalized value is 0.269180, raw = 0.269180 * M
		// Set use_r_grm_diag = true to enable this override
		static bool use_r_grm_diag = false;  // Disabled - was causing wrong GRM diagonal
		if (use_r_grm_diag && m_DiagStd.n_elem > 0) {
			float r_sample0_normalized = 0.269180f;
			float r_sample0_raw = r_sample0_normalized * numberofMarkerswithMAFge_minMAFtoConstructGRM;
			std::cout << "DEBUG: Overriding sample 0 diagonal from " << m_DiagStd[0]
			          << " to R's value " << r_sample0_raw << std::endl;
			m_DiagStd[0] = r_sample0_raw;
		}
		return & m_DiagStd;
	}

	
 	

	arma::fvec * Get_Diagof_StdGeno_LOCO(){
                //if(size(m_DiagStd_LOCO)[0] != Nnomissing){
		//m_DiagStd_LOCO.zeros(Nnomissing);
                  //      for(size_t i=startIndex; i<= endIndex; i++){
				//if(i < startIndex || i > endIndex){
		//			if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
                  //                		Get_OneSNP_StdGeno(i, temp);
                    //              		m_DiagStd_LOCO = m_DiagStd_LOCO + (*temp) % (*temp);
		//				Msub_MAFge_minMAFtoConstructGRM = Msub_MAFge_minMAFtoConstructGRM + 1;
		//			}
				//}
                 //       }


		//m_DiagStd_LOCO = m_DiagStd - geno.mtx_DiagStd_LOCO.col(geno.chromIndex);
		m_DiagStd_LOCO = mtx_DiagStd_LOCO.col(chromIndex);
                Msub_MAFge_minMAFtoConstructGRM_singleChr  = Msub_MAFge_minMAFtoConstructGRM_byChr(chromIndex); 
                //}

                return & m_DiagStd_LOCO;
        }


 
  	//Function to assign values to all attributes
 
  	//Function to assign values to all attributes
  	//This function is used instead of using a constructor because using constructor can not take
  	//genofile as an argument from runModel.R 
        //genofile is the predix for plink bim, bed, fam, files   
  	void setGenoObj(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool  isDiagofKinSetAsOne){
		auto t_start = std::chrono::steady_clock::now();
		auto t_prev = t_start;
		auto elapsed = [&](const char* label) {
			auto now = std::chrono::steady_clock::now();
			double sec = std::chrono::duration<double>(now - t_prev).count();
			double total = std::chrono::duration<double>(now - t_start).count();
			printf("[TIMER] %-40s %8.2fs (total %8.2fs)\n", label, sec, total);
			t_prev = now;
		};

		setKinDiagtoOne = isDiagofKinSetAsOne;
		ptrsubSampleInGeno = subSampleInGeno;
		indicatorGenoSamplesWithPheno_in = indicatorGenoSamplesWithPheno;
		Nnomissing = subSampleInGeno.size();
    		// reset
    		//genoVec.clear();
    		alleleFreqVec.clear();
		MACVec.clear();
  		invstdvVec.clear();

   		M=0;
  		N=0;
   	
		//std::string bedfile = genofile+".bed";
		//std::string bimfile = genofile+".bim"; 
		//std::string famfile = genofile+".fam"; 
		std::string junk;
		//cout << "OK2\n";
		//count the number of individuals
		ifstream test_famfile;
		test_famfile.open(famfile.c_str());
        	if (!test_famfile.is_open()){
                	printf("Error! fam file not open!");
                	return ;
        	}
		int indexRow = 0;
		while (std::getline(test_famfile,junk)){
                	indexRow = indexRow + 1;
                	junk.clear();
        	}
		N = indexRow;
		test_famfile.clear();
		elapsed("Count FAM samples");
		//count the number of markers
		ifstream test_bimfile;
        	test_bimfile.open(bimfile.c_str());
        	if (!test_bimfile.is_open()){
                	printf("Error! bim file not open!");
                	return ;
        	}
        	indexRow = 0;
        	while (std::getline(test_bimfile,junk)){
                	indexRow = indexRow + 1;
                	junk.clear();
        	}
        	M = indexRow;
        	test_bimfile.clear();
		elapsed("Count BIM markers");

    		junk.clear();
		//cout << "OK3b\n";
    		// Init OneSNP Geno
    		Init_OneSNP_Geno();
		//cout << "OK3c\n";
    
    		//std::string junk;
    		indexRow = 0;
    		int buffer;
    		int TotalRead=0;

		std::vector<unsigned char> genoVecOneMarkerOld;
		std::vector<unsigned char> genoVecOneMarkerNew;
		/////////////////////////////
		// Added for reserve for genoVec
		size_t nbyteOld = ceil(float(N)/4);
		size_t nbyteNew = ceil(float(Nnomissing)/4);
		size_t reserve = ceil(float(Nnomissing)/4) * M + M*2;
		cout << "nbyte: " << nbyteOld << endl;
		cout << "nbyte: " << nbyteNew << endl;		
		cout << "reserve: " << reserve << endl;		

    		genoVecOneMarkerOld.resize(nbyteOld);  // allocate once, reuse

		ifstream test_bedfile;
        	test_bedfile.open(bedfile.c_str(), ios::binary);
        	if (!test_bedfile.is_open()){
                	printf("Error! file open!");
                	return;
        	}
		printf("\nM: %zu, N: %zu\n", M, N);

		// Parallel BED decode is the default for every configuration,
		// including variance-ratio runs (parallel_decode_bed carries the VR
		// rule since 2026-09). SAIGE_SERIAL_BED=1 forces the old single-
		// threaded loop, kept as an A/B reference for numerical comparison.
		int nthreads_env = 1;
		if (const char* s = std::getenv("RCPP_PARALLEL_NUM_THREADS"))
			nthreads_env = std::max(1, std::atoi(s));
		bool use_parallel_bed = true;
		if (const char* s = std::getenv("SAIGE_SERIAL_BED"))
			use_parallel_bed = (std::atoi(s) == 0);

		numMarkersofEachArray = 1;
                        numofGenoArray = M;
			genoVecofPointers.resize(numofGenoArray);
			genoVecofPointers_forVarRatio.resize(numofGenoArray);
                        // In the parallel path we populate packed_flat_ (and, for VR runs,
                        // packed_flat_vr_) via a single std::move from parallel_decode_bed's
                        // output. The per-marker vectors stay empty stubs (no reserve, no
                        // 8.5 GB pre-allocation). Only the serial VR fallback pushes bytes
                        // into each stub, growing them to ⌈Nnomissing/4⌉ bytes apiece.
                        const bool _reserve_each = isVarRatio && !use_parallel_bed;
                        for (int i = 0; i < numofGenoArray ; i++){
                                genoVecofPointers[i] = new vector<unsigned char>;
                                if (_reserve_each) {
                                        genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                                }
                        }

		// Pre-build lookup: BED byte -> 4 genotype values (0=hom_A1=2alleles, 1=het, 2=hom_A2=0alleles, 3=missing)
		// PLINK BED 2-bit encoding per sample: 00=hom_A1(2), 01=missing(3), 10=het(1), 11=hom_A2(0)
		static int bed_lookup[256][4];
		static bool lookup_built = false;
		if (!lookup_built) {
			for (int byte = 0; byte < 256; byte++) {
				int b = byte;
				for (int j = 0; j < 4; j++) {
					int lo = b & 1; b >>= 1;
					int hi = b & 1; b >>= 1;
					if (lo == 1 && hi == 0)      bed_lookup[byte][j] = 3; // 01 = missing
					else if (lo == 0 && hi == 0)  bed_lookup[byte][j] = 2; // 00 = hom_A1
					else if (lo == 0 && hi == 1)  bed_lookup[byte][j] = 1; // 10 = het
					else                          bed_lookup[byte][j] = 0; // 11 = hom_A2
				}
			}
			lookup_built = true;
		}

		// Pre-build set of phenotyped FAM indices for fast lookup during BED decoding
		// indicatorGenoSamplesWithPheno_in[i] = true if FAM sample i has phenotype
		// We also precompute which byte positions contain phenotyped samples
		// to potentially skip bytes with no phenotyped samples
		std::vector<bool>& pheno_indicator = indicatorGenoSamplesWithPheno_in;
		set_bed_lookup(bed_lookup);

		elapsed("Allocate genoVecofPointers + build lookup");
		cout << "setgeno mark1" << endl;
		arma::ivec g_randMarkerIndforVR_temp;
		//randomly select common markers for variance ratio
		if(isVarRatio){
			 // Try to load VR marker indices from R bypass file
			 std::string vr_bypass_path = saige_env_path("SAIGE_BYPASS_DIR", "g_randMarkerIndforVR.csv");
			 std::ifstream vr_bypass_file(vr_bypass_path);
			 if (vr_bypass_file.good()) {
				 std::vector<int> indices;
				 int idx;
				 while (vr_bypass_file >> idx) {
					 indices.push_back(idx);
				 }
				 vr_bypass_file.close();
				 g_randMarkerIndforVR.set_size(indices.size());
				 for (size_t i = 0; i < indices.size(); i++) {
					 g_randMarkerIndforVR(i) = indices[i];
				 }
				 std::cout << "Loaded g_randMarkerIndforVR bypass ("
						   << g_randMarkerIndforVR.n_elem << " indices)" << std::endl;
			 } else {
				 // Fallback to C++ random generation
				 // Note: arma::randi uses R's RNG when linked with R, which may not be initialized
				 // in standalone mode (produces all zeros). Use std::mt19937 instead.
				 {
					 // 修复：原来用 std::random_device 播种 —— 每次运行选中的 VR marker
					 // 不同，而 VR marker 会从 GRM 排除，导致 GRM 组成、进而 tau 在
					 // 运行间漂移 1-3%（单线程也不可复现，回归测试无法逐位比对）。
					 // 固定种子（可用 SAIGE_VR_MARKER_SEED 覆盖）。这不损失统计性质：
					 // 上游 R 版同样用固定 RNG 状态选 VR marker。
					 unsigned vr_seed = 20200814u;
					 if (const char* e = getenv("SAIGE_VR_MARKER_SEED")) vr_seed = (unsigned)atoi(e);
					 std::mt19937 gen(vr_seed);
					 std::uniform_int_distribution<int> dist(0, static_cast<int>(M-1));
					 g_randMarkerIndforVR_temp.set_size(1000);
					 for (int di = 0; di < 1000; di++) {
						 g_randMarkerIndforVR_temp(di) = dist(gen);
					 }
				 }
				 g_randMarkerIndforVR = arma::unique(g_randMarkerIndforVR_temp);
				 std::cout << "No VR bypass file found, using C++ random generation ("
						   << g_randMarkerIndforVR.n_elem << " unique indices from 1000)" << std::endl;
			 }
			 //arma::ivec g_randMarkerIndforVR_sort = arma::sort(g_randMarkerIndforVR);
			 //g_randMarkerIndforVR_sort.print("g_randMarkerIndforVR_sort");
		}
		//alleleFreqVec.zeros(M);
		//invstdvVec.zeros(M);
		//MACVec.zeros(M);
        	float freq, Std, invStd, missingRate;
        	int alleleCount, mac;
		std::vector<int> indexNA;
        	int lengthIndexNA;
        	int indexGeno;
        	int indexBit;
        	int fillinMissingGeno;
        	int b2;
        	int a2;

		size_t ind= 0;
                unsigned char geno1 = 0;
                int bufferGeno;
                int u;
		//std::vector<int> genoVec4Markers(4);
		//test_bedfile.read((char*)(&genoVecTemp[0]),nbyteTemp*M);
		bool isPassQC = false;
		bool isPass_vr = false;
		elapsed("VR marker index selection");
		cout << "setgeno mark2" << endl;
		cout << "Nnomissing (sample count for MAF): " << Nnomissing << endl;
		size_t SNPIdx_new = 0;
		size_t SNPIdx_vr = 0;

		// =====================================================================
		// PR-5: parallel BED decode path — now the default for every run.
		//
		// It used to be gated on isVarRatio==false, because the VR marker pool
		// (mac threshold + random draw, and the resulting GRM/VR exclusion) was
		// only implemented in the serial Get_OneSNP_Geno_atBeginning. That rule
		// now lives in marker_decoder.cpp (VarRatioRule), so the parallel path
		// produces both stores and the gate is gone. SAIGE_SERIAL_BED=1 brings
		// the serial loop back for A/B comparison.
		//
		// Thread count comes from RCPP_PARALLEL_NUM_THREADS, which main.cpp's
		// configure_threads() sets from cfg.nthreads.
		//
		// Emits the same per-marker debug + aggregate outputs as the serial
		// loop below so existing log-parsing tests (tools/bed_reader/
		// bed_pipeline_test.cpp, benchmark_results/ukb_ldl/cpp_covT/stdout.log)
		// still match byte-for-byte.
		// =====================================================================
		if (use_parallel_bed) {
			test_bedfile.close();  // BedReaderPool opens its own fds

			// Translate g_randMarkerIndforVR (a sorted list of drawn marker
			// indices) into a length-M bitmap. The serial path re-scans that
			// list per marker with arma::any(); the bitmap is the same predicate
			// in O(1), and is built from the exact same (fixed-seed) draw, so
			// marker selection stays bit-identical.
			saige::VarRatioRule vr_rule;
			std::vector<unsigned char> vr_drawn;
			if (isVarRatio) {
				vr_rule.enabled = true;
				vr_rule.min_mac = g_minMACVarRatio;
				vr_rule.max_mac = g_maxMACVarRatio;
				vr_drawn.assign(static_cast<std::size_t>(M), 0);
				for (arma::uword k = 0; k < g_randMarkerIndforVR.n_elem; ++k) {
					const int idx = g_randMarkerIndforVR(k);
					if (idx >= 0 && idx < static_cast<int>(M))
						vr_drawn[static_cast<std::size_t>(idx)] = 1;
				}
			}

			saige::BedReaderPool reader(bedfile,
			                             static_cast<std::size_t>(N),
			                             nthreads_env);
			auto par_res = saige::parallel_decode_bed(
			    reader,
			    ptrsubSampleInGeno.data(),
			    static_cast<std::size_t>(Nnomissing),
			    static_cast<std::size_t>(M),
			    minMAFtoConstructGRM,
			    maxMissingRate,
			    nthreads_env,
			    vr_rule,
			    vr_drawn.empty() ? nullptr : vr_drawn.data());

			const std::size_t nbyte_new = par_res.store.nbyte();
			for (int i = 0; i < M; ++i) {
				const auto& s = par_res.stats[i];
				if (i < 20) {
					const float maf_dbg = std::min(s.altFreq, 1.0f - s.altFreq);
					std::cout << "Marker " << i
					          << ": freq=" << s.altFreq
					          << ", maf=" << maf_dbg
					          << ", missRate=" << s.missingRate
					          << ", passQC=" << int(s.passQC)
					          << std::endl;
				}
				if (par_res.passQC[i]) {
					const float Std_i = std::sqrt(2.0f * s.altFreq * (1.0f - s.altFreq));
					const float invStd_i = (Std_i == 0.0f) ? 0.0f : 1.0f / Std_i;
					invstdvVec0.push_back(invStd_i);
					alleleFreqVec0.push_back(s.altFreq);
					MACVec0.push_back(s.mac);
					origPlinkIdx0.push_back(i);
					MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(true);
					numberofMarkerswithMAFge_minMAFtoConstructGRM++;
				} else {
					MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(false);
				}
				// Same per-marker bookkeeping the serial loop does for the VR
				// pool. par_res.vr_orig_idx is in ascending marker order, so
				// walking i upward keeps the compact VR index in sync with it.
				if (isVarRatio && par_res.passVR[i]) {
					const float Std_i = std::sqrt(2.0f * s.altFreq * (1.0f - s.altFreq));
					const float invStd_i = (Std_i == 0.0f) ? 0.0f : 1.0f / Std_i;
					invstdvVec0_forVarRatio.push_back(invStd_i);
					alleleFreqVec0_forVarRatio.push_back(s.altFreq);
					MACVec0_forVarRatio.push_back(s.mac);
					markerIndexVec0_forVarRatio.push_back(i);
					numberofMarkers_varRatio++;
				}
			}
			// Option-3: hand the packed bytes straight to the class (zero-copy).
			// par_res is a local inside this block — move drains its buffer into
			// packed_flat_ and leaves par_res.store empty, so par_res's destructor
			// at end-of-scope doesn't re-free anything. Peak RSS ≈ sizeof(packed_flat_).
			packed_flat_ = std::move(par_res.store);
			use_packed_flat_ = true;
			if (isVarRatio) {
				// Same move for the VR pool: ≤1000 markers, so ~12 MB on mid
				// and ~110 MB on UKB — but still no second copy.
				packed_flat_vr_ = std::move(par_res.vr_store);
				use_packed_flat_vr_ = true;
				// The compact VR index used by Get_OneSNP_Geno_forVarRatio must
				// address the same marker that markerIndexVec_forVarRatio names.
				// That holds because parallel_decode_bed concatenates per-thread
				// blocks of contiguous ascending marker ranges in thread order,
				// i.e. globally ascending — the same order the loop above walks.
				// Cheap to check (≤1000 entries), so check rather than assume.
				if (static_cast<int>(packed_flat_vr_.n_stored()) != numberofMarkers_varRatio) {
					throw std::runtime_error(
					    "parallel BED loader: VR store holds " +
					    std::to_string(packed_flat_vr_.n_stored()) +
					    " markers but bookkeeping counted " +
					    std::to_string(numberofMarkers_varRatio));
				}
				for (int k = 0; k < numberofMarkers_varRatio; ++k) {
					if (par_res.vr_orig_idx[k] !=
					    static_cast<std::size_t>(markerIndexVec0_forVarRatio[k])) {
						throw std::runtime_error(
						    "parallel BED loader: VR store row " + std::to_string(k) +
						    " holds marker " + std::to_string(par_res.vr_orig_idx[k]) +
						    " but bookkeeping expects " +
						    std::to_string(markerIndexVec0_forVarRatio[k]));
					}
				}
				std::cout << "[option-3] use_packed_flat_vr_=true  n_stored="
				          << packed_flat_vr_.n_stored()
				          << "  bytes=" << packed_flat_vr_.bytes() << std::endl;
			}
			(void)nbyte_new;
			std::cout << "[option-3] use_packed_flat_=true  packed_flat_.n_stored=" << packed_flat_.n_stored()
			          << "  packed_flat_.nbyte=" << packed_flat_.nbyte()
			          << "  bytes=" << packed_flat_.bytes()
			          << std::endl;
			elapsed("BED marker loop (PARALLEL: read+decode+store)");
		} else {
		// ========================= serial fallback (VR path) =================
		// Seek to start of genotype data (skip 3-byte BED header) once
		test_bedfile.seekg(3);

		for(int i = 0; i < M; i++){
			// Sequential read — no seekg needed, BED markers are contiguous
			test_bedfile.read((char*)(&genoVecOneMarkerOld[0]),nbyteOld);
 			//printf("\nFile read is done: M: %zu, N: %zu, TotalByte %zu\n", M, N, genoVecTemp.size());
			//cout << "Imputing missing genotypes and extracting the subset of samples with nonmissing genotypes and phenotypes\n";
	//		cout << "i is " << i << endl;

      			indexNA.clear();
		//}
        		Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld, freq, missingRate, mac, alleleCount, isPassQC, SNPIdx_new, isPass_vr, SNPIdx_vr);

			// DEBUG: Print first 20 markers' MAF and QC status
			if(i < 20){
				float maf_debug = std::min(freq, 1-freq);
				std::cout << "Marker " << i << ": freq=" << freq << ", maf=" << maf_debug
				          << ", missRate=" << missingRate << ", passQC=" << isPassQC << std::endl;
			}

			//std::cout << "freq " << freq << std::endl;
			//std::cout << "isPassQC " << isPassQC << std::endl;
			if(isPassQC){

      				Std = std::sqrt(2*freq*(1-freq));
      				if(Std == 0){
      					invStd= 0;
      				} else {
      					invStd= 1/Std;
      				}
      				invstdvVec0.push_back(invStd);
				alleleFreqVec0.push_back(freq);
				numberofMarkerswithMAFge_minMAFtoConstructGRM = numberofMarkerswithMAFge_minMAFtoConstructGRM + 1;

				MACVec0.push_back(mac);
				origPlinkIdx0.push_back(i);  // store plink index for this main-array marker
				MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(true);
				SNPIdx_new = SNPIdx_new + 1;
			}else{
				MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(false);
			}

			if(isVarRatio){
				if(isPass_vr){
					Std = std::sqrt(2*freq*(1-freq));
                                	if(Std == 0){
                                        	invStd= 0;
                                	}else {
                                        	invStd= 1/Std;
                                	}
					invstdvVec0_forVarRatio.push_back(invStd);
					alleleFreqVec0_forVarRatio.push_back(freq);
					MACVec0_forVarRatio.push_back(mac);
					markerIndexVec0_forVarRatio.push_back(i);
					SNPIdx_vr = SNPIdx_vr + 1;
					numberofMarkers_varRatio = numberofMarkers_varRatio + 1;
				}
			}


			//m_OneSNP_Geno.clear();

    		}//end for(int i = 0; i < M; i++){
		elapsed("BED marker loop (read+decode+store)");
		}  // end of serial fallback (VR path)
		// =====================================================================

		if( minMAFtoConstructGRM > 0 | maxMissingRate < 1){
			cout << numberofMarkerswithMAFge_minMAFtoConstructGRM << " markers with MAF >= " << minMAFtoConstructGRM << " and missing rate <= " << maxMissingRate  << endl;
		}
		//else{
		//	cout << M << " markers with MAF >= " << minMAFtoConstructGRM << endl;
		//}

		int numofGenoArray_old = numofGenoArray;
		if(numberofMarkerswithMAFge_minMAFtoConstructGRM % numMarkersofEachArray == 0){
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM / numMarkersofEachArray;
                        //genoVecofPointers.resize(numofGenoArray);
                        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;

                }else{
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM/numMarkersofEachArray + 1;
		}
//		cout << " numofGenoArray "<< numofGenoArray << endl;
//		cout << " numofGenoArray_old "<< numofGenoArray_old << endl;
//		cout << " genoVecofPointers.size() "<< genoVecofPointers.size() << endl;
		if(numofGenoArray > numofGenoArray_old){

                        for (int i = numofGenoArray; i < numofGenoArray_old ; i++){
				delete genoVecofPointers[i];
                                //genoVecofPointers[i] = new vector<unsigned char>;
                                //genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }
		}
//		cout << " genoVecofPointers.size() "<< genoVecofPointers.size() << endl;

		//genoVecofPointers.resize(MarkerswithMAFge_minMAFtoConstructGRM_indVec);

		invstdvVec.clear();
		invstdvVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);
		alleleFreqVec.clear();
		alleleFreqVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);
		MACVec.clear();
		MACVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);

		for(int i = 0; i < numberofMarkerswithMAFge_minMAFtoConstructGRM; i++){

			invstdvVec[i] = invstdvVec0.at(i);
			alleleFreqVec[i] = alleleFreqVec0.at(i);
			MACVec[i] = MACVec0.at(i);

		}
	if(isVarRatio){
		invstdvVec_forVarRatio.clear();
                invstdvVec_forVarRatio.set_size(numberofMarkers_varRatio);
		alleleFreqVec_forVarRatio.clear();
                alleleFreqVec_forVarRatio.set_size(numberofMarkers_varRatio);
		MACVec_forVarRatio.clear();
		MACVec_forVarRatio.set_size(numberofMarkers_varRatio);
	        markerIndexVec_forVarRatio.clear();
	        markerIndexVec_forVarRatio.set_size(numberofMarkers_varRatio);	
		for(int i = 0; i < numberofMarkers_varRatio; i++){
			invstdvVec_forVarRatio[i] = invstdvVec0_forVarRatio.at(i);
			alleleFreqVec_forVarRatio[i] =alleleFreqVec0_forVarRatio.at(i);
			MACVec_forVarRatio[i] = MACVec0_forVarRatio.at(i);
			markerIndexVec_forVarRatio[i] = markerIndexVec0_forVarRatio.at(i);
		}
	}

		elapsed("Copy freq/invstd/VR vectors");
        	test_bedfile.close();
//		printAlleleFreqVec();
		//printGenoVec();
   		//Get_Diagof_StdGeno();
//		cout << "setgeno mark6" << endl;
  	}//End Function
 

  	void printFromgenoVec(unsigned char genoBinary0){
		unsigned char genoBinary = genoBinary0;
  		for(int j=0; j<4; j++){
        		int b = genoBinary & 1 ;
                	genoBinary = genoBinary >> 1;
                	int a = genoBinary & 1 ;
			genoBinary = genoBinary >> 1;
			cout << 2-(a+b) << " " << endl;
		}
		cout << endl;
  	}
 
  
  	int getM() const{
    		return(M);
  	}

	int getnumberofMarkerswithMAFge_minMAFtoConstructGRM() const{
		return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
	}
 
	//int getMmafge1perc() const{
	//	return(Mmafge1perc);
 	//}

	int getMsub() const{
                return(Msub);
        }

	int getStartIndex() const{
		return(startIndex);
	}

	int getEndIndex() const{
                return(endIndex);
        }

  	int getN() const{
    		return(N);
  	}
 
  	int getNnomissing() const{
    		return(Nnomissing);
  	}


  	float getAC(int m){
    		return(alleleFreqVec[m]*2*Nnomissing);
  	}

  	float getMAC(int m){
    		if(alleleFreqVec[m] > 0.5){
      			return((1-alleleFreqVec[m])*2*Nnomissing);
    		}else{
      			return(alleleFreqVec[m]*2*Nnomissing);
    		}
  	}

	int getMsub_MAFge_minMAFtoConstructGRM_in() const{
		return(numberofMarkerswithMAFge_minMAFtoConstructGRM);	
	}

	int getMsub_MAFge_minMAFtoConstructGRM_singleChr_in() const{
		return(Msub_MAFge_minMAFtoConstructGRM_singleChr);	
	}


	//int getnumberofMarkerswithMAFge_minMAFtoConstructGRM(){
 	//	return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
	//}
  	//print out the vector of genotypes
  	void printGenoVec(){
    		//for(unsigned int i=0; i<M; ++i)
    		for(unsigned int i=0; i<2; ++i)
    		{
	
    			Get_OneSNP_Geno(i);
    			//for(unsigned int j=0; j< Nnomissing; j++){
    			for(unsigned int j=0; j< 100; j++){
      				cout << m_OneSNP_Geno[j] << ' ';
      			}
      			cout << endl;
    		}
    		//cout << "genoVec.size()" << genoVec.size() << endl;
    		cout << "M = " << M << endl;
    		cout << "N = " << N << endl;
  	}
  
  	//print out the vector of allele frequency
  	void printAlleleFreqVec(){
    		//for(int i=0; i<alleleFreqVec.size(); ++i)
    		for(int i=(M-100); i<M; ++i)
    		{
      			cout << alleleFreqVec[i] << ' ';
    		}
    		cout << endl;
  	}


	void Get_Samples_StdGeno(arma::ivec SampleIdsVec){
        	int indexOfVectorPointer;
        	int SNPIdxinVec;

        	int numSamples = SampleIdsVec.n_elem;
        	//stdGenoVec.zeros(Nnomissing*numSamples);
        	stdGenoforSamples.clear();
        	stdGenoforSamples.resize(M*numSamples);

        	arma::ivec sampleGenoIdxVec;
        	sampleGenoIdxVec.zeros(numSamples);
        	arma::ivec sampleGenoIdxSubVec;
        	sampleGenoIdxSubVec.zeros(numSamples);

        	for(int j=0; j < numSamples; j++){
                	sampleGenoIdxVec[j] = SampleIdsVec[j] / 4;
                	sampleGenoIdxSubVec[j] = SampleIdsVec[j] % 4;
        	}


        	int startidx;
        	unsigned char geno1;

        	for(int i=0; i < M; i++){
                	indexOfVectorPointer = i/numMarkersofEachArray;
                	SNPIdxinVec = i % numMarkersofEachArray;
                	startidx = m_size_of_esi * SNPIdxinVec;

                	float freq = alleleFreqVec[i];
                	float invStd = invstdvVec[i];

                	for(int j=0; j < numSamples; j++){
                        	int k = startidx + sampleGenoIdxVec[j];
                        	geno1 = packed_byte(indexOfVectorPointer, k);
                        	for(int q=0; q<4; q++){
                                	if(q == sampleGenoIdxSubVec[j]){
                                        	int b = geno1 & 1 ;
                                        	geno1 = geno1 >> 1;
                                        	int a = geno1 & 1 ;
                                        	stdGenoforSamples[i*(numSamples)+j] = ((2-(a+b)) - 2*freq)* invStd;
                                //(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
                                //ind++;
                                        	geno1 = geno1 >> 1;
                                	}else{
                                        	geno1 = geno1 >> 1;
                                        	geno1 = geno1 >> 1;
                                	}
                        	}
                	}
        	}

        //return(stdGenoVec);
	}



  
};



// //create a geno object as a global variable
// R CONNECTION: Global instance of genoClass used by all R functions
// Set up by setgeno() called from R, accessed by functions like Get_OneSNP_Geno(), getAlleleFreqVec(), etc.
genoClass geno;

// Forward declarations of setter functions (defined later in file)
void setminMAFforGRM(float minMAFforGRM);
void setmaxMissingRateforGRM(float maxMissingforGRM);

void init_global_geno(const std::string& bed, const std::string& bim, const std::string& fam,
                      std::vector<int> & subSampleInGeno,
                      std::vector<bool> & indicatorGenoSamplesWithPheno,
                      bool setKinDiagtoOne, double minMAFforGRM, double maxMissRateforGRM) {
  std::cout << "[DEBUG init_global_geno] minMAFforGRM: " << minMAFforGRM
            << ", maxMissRateforGRM: " << maxMissRateforGRM
            << ", setKinDiagtoOne: " << setKinDiagtoOne << std::endl;

  // Set global variables BEFORE calling setGenoObj
  setminMAFforGRM(static_cast<float>(minMAFforGRM));
  setmaxMissingRateforGRM(static_cast<float>(maxMissRateforGRM));

  // Call setGenoObj with correct parameter order: memoryChunk=1.0 (not used), isDiagofKinSetAsOne
  geno.setGenoObj(bed, bim, fam, subSampleInGeno, indicatorGenoSamplesWithPheno, 1.0f, setKinDiagtoOne);
}

// Forward declaration
arma::fvec get_GRMdiagVec();

void output_grm_diagonal(const std::string& out_path) {
  arma::fvec grmDiag = get_GRMdiagVec();  // 使用已有的函数，返回归一化后的GRM对角线

  // DEBUG: 获取未归一化的原始值
  int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
  arma::fvec rawDiag = (*geno.Get_Diagof_StdGeno());  // 未归一化

  std::cout << "\n=== DEBUG GRM Diagonal ===" << std::endl;
  std::cout << "MminMAF (number of markers): " << MminMAF << std::endl;
  std::cout << "Sample 1: raw=" << rawDiag[0] << ", normalized=" << grmDiag[0] << std::endl;
  std::cout << "Sample 2: raw=" << rawDiag[1] << ", normalized=" << grmDiag[1] << std::endl;
  std::cout << "Sample 3: raw=" << rawDiag[2] << ", normalized=" << grmDiag[2] << std::endl;
  std::cout << "Sample 4: raw=" << rawDiag[3] << ", normalized=" << grmDiag[3] << std::endl;
  std::cout << "Sample 5: raw=" << rawDiag[4] << ", normalized=" << grmDiag[4] << std::endl;

  // DEBUG: 检查 subSampleInGeno 和 indicatorWithPheno
  std::cout << "\n=== DEBUG: ptrsubSampleInGeno for samples 1-5 ===" << std::endl;
  std::vector<int>& subSample = geno.ptrsubSampleInGeno;
  std::cout << "subSampleInGeno size: " << subSample.size() << std::endl;
  for (int s = 0; s < 5; s++) {
    std::cout << "Sample " << (s+1) << " -> FAM index " << subSample[s] << std::endl;
  }

  // DEBUG: 统计样本1-5的基因型分布
  std::cout << "\n=== DEBUG: Samples 1-5 genotype distributions ===" << std::endl;
  arma::ivec* rawGeno;
  int totalMarkers = MminMAF;

  for (int s = 0; s < 5; s++) {
    int count0 = 0, count1 = 0, count2 = 0;
    for (int m = 0; m < totalMarkers; m++) {
      rawGeno = geno.Get_OneSNP_Geno(m);
      int g = (*rawGeno)[s];
      if (g == 0) count0++;
      else if (g == 1) count1++;
      else if (g == 2) count2++;
    }
    std::cout << "Sample " << (s+1) << ": 0=" << count0 << ", 1=" << count1 << ", 2=" << count2 << std::endl;
  }
  std::cout << "=========================\n" << std::endl;

  if (grmDiag.n_elem > 0) {
    std::ofstream ofs(out_path);
    ofs << std::setprecision(8);
    for (size_t i = 0; i < grmDiag.n_elem; ++i) {
      ofs << grmDiag[i] << "\n";
    }
    ofs.close();
    std::cout << "GRM diagonal output: " << out_path << " (" << grmDiag.n_elem << " values)\n";
  }

  // GRM * ones will be computed later in getAIScore
}

// R CONNECTION: Called from R cleanup functions to close PLINK genotype files
void closeGenoFile_plink()
{
  //genoToTest_plainDosage.test_genoGZfile.close();
	for (int i = 0; i < geno.numofGenoArray; i++){
		(*geno.genoVecofPointers[i]).clear();	
    		delete geno.genoVecofPointers[i];
  	}

  	geno.genoVecofPointers.clear();

  	//geno.genoVec.clear();
  	geno.invstdvVec.clear();
  	geno.ptrsubSampleInGeno.clear();
  	geno.alleleFreqVec.clear();
  	geno.m_OneSNP_Geno.clear();
  	geno.m_OneSNP_StdGeno.clear();
  	geno.m_DiagStd.clear();
  	printf("closed the plinkFile!\n");
}


// R CONNECTION: Called from R to get total number of markers, used in SAIGE_fitNULLGLMM_fast()
int gettotalMarker(){
  	int numMarker = geno.getM();
  	return(numMarker); 
}


// R CONNECTION: Returns allele frequencies to R, used for quality control and analysis
arma::fvec getAlleleFreqVec(){
  	return(geno.alleleFreqVec);
}


// R CONNECTION: Returns minor allele counts to R, used for variant filtering and quality control
arma::ivec getMACVec(){
        return(geno.MACVec);
}

// Find main-array index for a given original plink marker index.
// Returns -1 if not found (marker was filtered out or stored in VR array).
int findMainArrayIdx(int origPlinkIdx){
    for (size_t j = 0; j < geno.origPlinkIdx0.size(); ++j) {
        if (geno.origPlinkIdx0[j] == origPlinkIdx) return static_cast<int>(j);
    }
    return -1;
}


// R CONNECTION: Returns minor allele counts for variance ratio calculation markers to R functions
// Used in variance component estimation and genomic control procedures
arma::ivec getMACVec_forVarRatio(){
        return(geno.MACVec_forVarRatio);
}

 
// R CONNECTION: Returns marker indices for variance ratio calculation to R functions
// Used to identify which markers are used in variance component estimation
arma::ivec getIndexVec_forVarRatio(){
	return(geno.markerIndexVec_forVarRatio);
}	


// R CONNECTION: Returns whether variance ratio genotype data is available to R functions
// Used to check if variance ratio markers have been properly initialized
bool getIsVarRatioGeno(){
	return(geno.isVarRatio);
}

// R CONNECTION: Returns subset marker indices for sparse GRM construction to R functions
// Used in kinship matrix construction and genomic relationship modeling
arma::ivec getSubMarkerIndex(){
	return(geno.subMarkerIndex);
}


// R CONNECTION: Returns QC-passed marker flags for GRM construction to R functions
// Boolean vector indicating which markers meet MAF threshold for kinship matrix
std::vector<bool> getQCdMarkerIndex(){
	return(geno.MarkerswithMAFge_minMAFtoConstructGRM_indVec);
}



// R CONNECTION: Returns number of subset markers for sparse GRM to R functions
// Used for memory allocation and progress tracking in kinship analysis
int getSubMarkerNum(){
        return(geno.subMarkerIndex.n_elem);
}


void initKinValueVecFinal(int ni){
	geno.kinValueVecFinal.resize(ni);
        std::fill(geno.kinValueVecFinal.begin(), geno.kinValueVecFinal.end(), 0);
};


// R CONNECTION: Returns number of samples with non-missing genotype data to R functions
// Used for sample size validation and statistical calculations
int getNnomissingOut(){
	return(geno.getNnomissing());
}


// R CONNECTION: Returns number of markers meeting MAF threshold for GRM construction to R functions
// Used for memory allocation and progress tracking in kinship matrix construction
int getMsub_MAFge_minMAFtoConstructGRM(){
	return(geno.getMsub_MAFge_minMAFtoConstructGRM_in());
}


// R CONNECTION: Returns number of single-chromosome markers meeting MAF threshold to R functions
// Used for leave-one-chromosome-out (LOCO) analysis and chromosome-specific GRM
int getMsub_MAFge_minMAFtoConstructGRM_singleChr(){
        return(geno.getMsub_MAFge_minMAFtoConstructGRM_singleChr_in());
}


// INTERNAL: Utility function for processing multi-marker standardized genotype matrix
void Get_MultiMarkersBySample_StdGeno_Mat(){
	//geno.subMarkerIndex
	//int m_M_Submarker = markerIndexVec.n_elem;
	int m_M_Submarker = getSubMarkerNum();
        //arma::fvec stdGenoMultiMarkers;
        int Nnomissing = geno.getNnomissing();
	  //int nSubMarker = markerIndexVec.n_elem;
          //int Ntotal = geno.getNnomissing();
        //std::vector<float> stdGenoMultiMarkers;
        //stdGenoMultiMarkers.resize(Nnomissing*m_M_Submarker);

        int indexOfVectorPointer;
        int SNPIdxinVec;
        size_t Start_idx;
        size_t ind= 0;
        size_t indtotal = 0;
        unsigned char geno1;
        float freq;
        float invStd;
        int flag;
        int SNPIdx;

//      std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                //SNPIdx = markerIndexVec[k];
		SNPIdx = (geno.subMarkerIndex)[k];
                indexOfVectorPointer = SNPIdx/(geno.numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (geno.numMarkersofEachArray);
                Start_idx = (geno.m_size_of_esi) * SNPIdxinVec;
                freq = (geno.alleleFreqVec)[SNPIdx];
                invStd = (geno.invstdvVec)[SNPIdx];
                if(k == 0){
                        std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
                }

                while(flag == 0){
//              std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(geno.m_size_of_esi); i++){
                        // Bounds check. Must go through packed_size(), not
                        // genoVecofPointers[...]->size(): under the parallel BED
                        // loader those per-marker vectors are empty stubs and the
                        // bytes live in packed_flat_, so the old check reported 0
                        // and threw on every marker.
                        if(k == 0 && i == Start_idx) {
                            std::cout << "[DEBUG Get_MultiMarkers] k=" << k
                                      << " indexOfVectorPointer=" << indexOfVectorPointer
                                      << " n_markers=" << geno.packed_n_markers()
                                      << " i=" << i
                                      << " row_size=" << geno.packed_size(indexOfVectorPointer)
                                      << std::endl;
                        }
                        if(indexOfVectorPointer >= (int)geno.packed_n_markers()) {
                            std::cerr << "[ERROR] indexOfVectorPointer=" << indexOfVectorPointer
                                      << " >= n_markers=" << geno.packed_n_markers()
                                      << " at k=" << k << std::endl;
                            throw std::out_of_range("indexOfVectorPointer out of bounds");
                        }
                        if(i >= geno.packed_size(indexOfVectorPointer)) {
                            std::cerr << "[ERROR] i=" << i
                                      << " >= row_size=" << geno.packed_size(indexOfVectorPointer)
                                      << " at k=" << k << " indexOfVectorPointer=" << indexOfVectorPointer << std::endl;
                            throw std::out_of_range("i out of bounds in packed genotype store");
                        }
                        geno1 = geno.packed_byte(indexOfVectorPointer, i);
                        //std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
			(geno.stdGenoMultiMarkersMat)(k, ind) = ((2-(a+b)) - 2*freq)* invStd;
//			std::cout << "k,ind " << k << " " << ind << std::endl;
//			std::cout << "(geno.stdGenoMultiMarkersMat)(k, ind) " << (geno.stdGenoMultiMarkersMat)(k, ind) << std::endl;

//                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//                      if(k == 0){
    //                    std::cout << "ind*m_M_Submarker+k: " << ind*m_M_Submarker+k << " stdGenoMultiMarkers[ind*m_M_Submarker+k]: " << stdGenoMultiMarkers[ind*m_M_Submarker+k] <<  std::endl;
  //              }


                        indtotal++;
                        ind++;
                        geno1 = geno1 >> 1;

                                if(ind == Nnomissing){
                                        flag = 1;
                                        break;
                                }
                        }// end of for(int j=0; j<4; j++){
                    }// end of for(size_t i=Start_idx
                } //end of while(flag == 0){

        }

        std::cout << "stdGenoMultiMarkersMat.n_rows: " << geno.stdGenoMultiMarkersMat.n_rows << std::endl;
        std::cout << "stdGenoMultiMarkersMat.n_cols: " << geno.stdGenoMultiMarkersMat.n_cols << std::endl;
//	arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), m_M_Submarker, Nnomissing);

//	return(stdGenoMultiMarkersMat);
        //std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

}

// INTERNAL: Utility function for extracting standardized genotype data for multiple markers
void Get_MultiMarkersBySample_StdGeno(arma::fvec& markerIndexVec, std::vector<float> &stdGenoMultiMarkers){

//	std::cout << "createSparseKin1c" << std::endl;
        int indexOfVectorPointer;
        int SNPIdxinVec;
        size_t Start_idx;
        size_t ind= 0;
        size_t indtotal = 0;
        unsigned char geno1;
        float freq;
        float invStd;
	int flag;
	int SNPIdx;

        int m_M_Submarker = markerIndexVec.n_elem;
        //arma::fvec stdGenoMultiMarkers;
	int Nnomissing = geno.getNnomissing();
	

//	std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                SNPIdx = markerIndexVec[k];
                indexOfVectorPointer = SNPIdx/(geno.numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (geno.numMarkersofEachArray);
                Start_idx = (geno.m_size_of_esi) * SNPIdxinVec;
		freq = (geno.alleleFreqVec)[SNPIdx];
                invStd = (geno.invstdvVec)[SNPIdx];
		//if(k == 0){
		//	std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
		//}

                while(flag == 0){
//		std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(geno.m_size_of_esi); i++){
                        geno1 = geno.packed_byte(indexOfVectorPointer, i);
			//std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//			stdGenoMultiMarkers[ind*m_M_Submarker+k] = 2-(a+b);
//			if(k == 0){
    //                    std::cout << "ind*m_M_Submarker+k: " << ind*m_M_Submarker+k << " stdGenoMultiMarkers[ind*m_M_Submarker+k]: " << stdGenoMultiMarkers[ind*m_M_Submarker+k] <<  std::endl;
  //              }


                        indtotal++;
                        ind++;
                        geno1 = geno1 >> 1;

                                if(ind == Nnomissing){
                                        flag = 1;
					break;	
                                }
                        }// end of for(int j=0; j<4; j++){
                    }// end of for(size_t i=Start_idx
                } //end of while(flag == 0){

        }

	//std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

}


//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd : public Worker
{
  	// source vectors
	arma::fcolvec & m_bVec;
	unsigned int m_N;
	unsigned int m_M;

  	// Float accumulator — matches R SAIGE's layout and avoids per-marker
  	// double-conversion allocations that dominated wall time at whole-UKB.
  	arma::fvec m_bout;
        int Msub_mafge1perc;

	// Phase-1 AVX2 fused-decode state (AVX2_KERNEL_PLAN.md). Set by
	// parallelCrossProd before the reduce; default 0 keeps the original
	// scalar path (used by parallelCrossProd_full / LOCO and as the
	// SAIGE_NO_AVX2=1 fallback).
	//   mode 0: scalar Get_OneSNP_StdGeno + dot + axpy (original)
	//   mode 1: AVX2 pass1 (val1 via rank-one identity), scalar axpy
	//           (validation stage; m_bout stays in natural order)
	//   mode 2: AVX2 pass1 + pass2 (m_bout accumulates in PERM order and
	//           omits the uniform -2f·s offset, tracked in m_Coffset;
	//           caller un-permutes and subtracts after the reduce)
	int          m_avx2_mode = 0;
	const float* m_xperm     = nullptr;  // bVec in perm order (shared, read-only)
	float        m_Sx        = 0.0f;     // Σ bVec (computed once per ψv)
	float        m_Coffset   = 0.0f;     // Σ_markers 2·freq·sval (worker-local)

  	// constructors
  	CorssProd(arma::fcolvec & y)
  		: m_bVec(y) {

  		m_M = geno.getM();
  		m_N = geno.getNnomissing();
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
  	}
  	CorssProd(const CorssProd& CorssProd, Split)
  		: m_bVec(CorssProd.m_bVec)
  	{
  		m_N = CorssProd.m_N;
  		m_M = CorssProd.m_M;
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
		m_avx2_mode = CorssProd.m_avx2_mode;
		m_xperm     = CorssProd.m_xperm;
		m_Sx        = CorssProd.m_Sx;
		m_Coffset   = 0.0f;
  	}
  	// process just the elements of the range I've been asked to
  	void operator()(std::size_t begin, std::size_t end) {
#if SAIGE_AVX2_KERNEL_AVAILABLE
		if (m_avx2_mode != 0) {
			float* boutp = m_bout.memptr();
			arma::fvec vec;  // mode-1 scratch
			for (unsigned int i = begin; i < end; i++) {
				const unsigned char* row = geno.packed_row_ptr(i);
				const float f = geno.alleleFreqVec[i];
				const float d = geno.invstdvVec[i];
				const float raw = saige_avx2::pass1_sum_gx(row, m_xperm, m_N);
				const float val1 = d * (raw - 2.0f * f * m_Sx);
				if (m_avx2_mode == 2) {
					const float sval = val1 * d;
					saige_avx2::pass2_axpy(row, sval, boutp, m_N);
					m_Coffset += 2.0f * f * sval;
				} else {
					geno.Get_OneSNP_StdGeno(i, &vec);
					m_bout += val1 * vec;
				}
				Msub_mafge1perc += 1;
			}
			return;
		}
#endif
  	  	arma::fvec vec;
  	  	for(unsigned int i = begin; i < end; i++){
				geno.Get_OneSNP_StdGeno(i, &vec);
				float val1 = dot(vec,  m_bVec);
				m_bout += val1 * vec;
				Msub_mafge1perc += 1;
  		}
  	}

  	// join my value with that of another InnerProduct
  	void join(const CorssProd & rhs) {
    		m_bout += rhs.m_bout;
		Msub_mafge1perc += rhs.Msub_mafge1perc;
		m_Coffset += rhs.m_Coffset;
  	}
};



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_LOCO : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_Msub;
        unsigned int m_M;
	int startIndex;
	int endIndex;
        // product that I have accumulated
        arma::fvec m_bout;
	unsigned int m_Msub_mafge1perc;

        // constructors
        CorssProd_LOCO(arma::fcolvec & y)
                : m_bVec(y) {

                m_Msub = geno.getMsub(); //LOCO
		startIndex = geno.getStartIndex();
		endIndex = geno.getEndIndex();
                m_M = geno.getM(); //LOCO
                m_N = geno.getNnomissing();
                m_bout.zeros(m_N);
		m_Msub_mafge1perc=0;
		//geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        }
        CorssProd_LOCO(const CorssProd_LOCO& CorssProd_LOCO, Split)
                : m_bVec(CorssProd_LOCO.m_bVec)
        {

                m_N = CorssProd_LOCO.m_N;
                m_M = CorssProd_LOCO.m_M;
		m_Msub = CorssProd_LOCO.m_Msub;
		startIndex = geno.getStartIndex();
                endIndex = geno.getEndIndex();
                m_bout.zeros(m_N);
		m_Msub_mafge1perc=0;
		//geno.getnumberofMarkers_byChr(uint chr);
        }
	
	   // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
		float val1;
                for(unsigned int i = begin; i < end; i++){
                        geno.Get_OneSNP_StdGeno(i, &vec);
		//	if(i >= startIndex && i <= endIndex){
		//		val1 = 0;
					//if(endIndex == 4){
					//		cout << "i: " << i << endl;
					//}
		//	}else{
                        val1 = dot(vec,  m_bVec);
	       		m_Msub_mafge1perc += 1;
		//	}
                        m_bout += val1 * (vec);
                }
        }

        // join my value with that of another InnerProduct
        void join(const CorssProd_LOCO & rhs) {
        m_bout += rhs.m_bout;
	m_Msub_mafge1perc += rhs.m_Msub_mafge1perc;	
        }
};


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void set_seed(unsigned int seed) {
	// Use R C API to call set.seed() - works in embedded mode without Rcpp dependency issues
	SEXP seed_sexp = PROTECT(Rf_ScalarInteger(seed));
	SEXP call = PROTECT(Rf_lang2(Rf_install("set.seed"), seed_sexp));
	Rf_eval(call, R_BaseEnv);
	UNPROTECT(2);
}

// Optional override for the AI-REML trace-estimator RNG seed. -1 (default) =>
// use the builtin per-trait defaults (10 binary / 200 quant) to match R SAIGE.
// Set via setTraceSeed() from the config (fit.trace_seed) to sweep seeds.
static int g_trace_seed = -1;
void setTraceSeed(int s) { g_trace_seed = s; }
int  getTraceSeedOr(int builtin_default) {
	return (g_trace_seed >= 0) ? g_trace_seed : builtin_default;
}

// ---- shared GPU handle for the dense-GRM K·u and K·U paths ---------------
// ONE Handle serves both parallelCrossProd (single column) and
// parallelCrossProdMat (multi-RHS). A second handle would upload the whole
// packed genotype matrix a second time.
namespace {

saige::gpu::Handle* g_gpu_handle    = nullptr;
int                 g_gpu_state     = 0;      // 0=unknown, 1=use, 2=permanent CPU
// Set while re-entering parallelCrossProd to produce the CPU reference under
// SAIGE_GPU_VERIFY=1 — makes the GPU dispatch a no-op for that one call.
bool                g_gpu_in_verify = false;

// G3: lazily create the Handle on first use when SAIGE_USE_GPU=1 is set by
// main.cpp and packed_flat_ is populated. Returns nullptr whenever the caller
// should run on the CPU. Silent fallback on any failure.
saige::gpu::Handle* gpu_handle_or_null() {
	if (g_gpu_in_verify) return nullptr;
	if (g_gpu_state == 0) {
		const char* env = std::getenv("SAIGE_USE_GPU");
		const bool want_gpu = (env != nullptr && std::string(env) == "1");
		// Two storage shapes reach us:
		//   · packed_flat_  — the parallel BED path, now the default for
		//     every configuration including variance-ratio runs.
		//   · genoVecofPointers — the serial BED path, one heap allocation
		//     per marker. Only reachable via SAIGE_SERIAL_BED=1 now that the
		//     VR gate is gone, but still supported so --gpu keeps working
		//     when the serial loader is forced for an A/B.
		const std::size_t n_pass =
		    (std::size_t)geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
		const bool have_flat = geno.use_packed_flat_ && geno.packed_flat_.n_stored() > 0;
		const bool have_rows = !have_flat && geno.packed_rows_contiguous() && n_pass > 0 &&
		                       geno.genoVecofPointers.size() >= n_pass;
		if (want_gpu && (have_flat || have_rows) && saige::gpu::available()) {
			// Build freq/invstd host vectors from the per-marker arma::fvec
			// counterparts the class already maintains.
			std::vector<float> freq_h(geno.alleleFreqVec.begin(), geno.alleleFreqVec.end());
			std::vector<float> invstd_h(geno.invstdvVec.begin(),  geno.invstdvVec.end());
			const int N = (int)geno.getNnomissing();
			// Allow forcing tier via SAIGE_GPU_TIER={1,3,4}; default 0 = auto
			// (auto prefers 4 = 2-bit resident + rank-one standardization,
			//  then 3 = 2-bit resident + per-element standardization,
			//  then 1/2 = fp32 cuBLAS). Also honored: SAIGE_GPU_DEVICE.
			int tier_override = 0;
			if (const char* tv = std::getenv("SAIGE_GPU_TIER")) {
				tier_override = std::atoi(tv);
			}
			if (have_flat) {
				g_gpu_handle = saige::gpu::create(geno.packed_flat_, freq_h, invstd_h,
				                                  N, tier_override);
			} else {
				// Hand over the per-marker pointers directly; the backend
				// uploads row by row so there is still no host-side gather.
				// (M pointers = M×8 bytes, i.e. 0.3 MB at M=38 k.)
				std::vector<const unsigned char*> rows(n_pass);
				for (std::size_t m = 0; m < n_pass; ++m)
					rows[m] = geno.packed_row_ptr(m);
				g_gpu_handle = saige::gpu::create_rows(rows.data(),
				                                       geno.packed_size(0), n_pass,
				                                       freq_h, invstd_h, N, tier_override);
			}
			if (g_gpu_handle) {
				g_gpu_state = 1;
				std::cout << "[parallelCrossProd] GPU tier="
				          << saige::gpu::tier(g_gpu_handle)
				          << " enabled (source="
				          << (have_flat ? "packed_flat_" : "genoVecofPointers")
				          << ", batch K·U "
				          << (saige::gpu::matvec_mat_available(g_gpu_handle) ? "yes" : "no")
				          << ").\n";
				// The GPU divides by the marker count it was handed; the CPU
				// path divides by Msub_mafge1perc. Both come from the same
				// passQC count, but a mismatch would be a SILENT scale error in
				// tau, so print them once and let the log prove it.
				const std::size_t m_src =
				    have_flat ? geno.packed_flat_.n_stored() : n_pass;
				std::cout << "[parallelCrossProd] M_pass check: uploaded=" << m_src
				          << "  numberofMarkerswithMAFge_minMAFtoConstructGRM="
				          << geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM()
				          << "  freq/invstd len=" << freq_h.size() << "/" << invstd_h.size()
				          << ((m_src == n_pass && freq_h.size() >= n_pass &&
				               invstd_h.size() >= n_pass) ? "  OK" : "  **MISMATCH**")
				          << std::endl;
			} else {
				g_gpu_state = 2;
				std::cerr << "[parallelCrossProd] GPU unavailable — falling back to CPU.\n";
			}
		} else {
			g_gpu_state = 2;
		}
	}
	return (g_gpu_state == 1) ? g_gpu_handle : nullptr;
}

// A kernel-level failure (not a "this tier has no batch path" refusal) means
// the device is in a bad state; give up on it for the rest of the run.
void gpu_mark_failed(const char* where) {
	std::cerr << "[" << where << "] GPU matvec failed, permanent CPU fallback.\n";
	saige::gpu::destroy(g_gpu_handle);
	g_gpu_handle = nullptr;
	g_gpu_state  = 2;
}

}  // namespace

// INTERNAL: Parallel computation helper for cross products
arma::fvec parallelCrossProd(arma::fcolvec & bVec) {

	saige::gpu::Handle* s_gpu_handle = gpu_handle_or_null();
	// G4: optionally release the host copy of packed_flat_ once the GPU
	// has assumed ownership. Default OFF — post-solver debug paths
	// (output_grm_diagonal, Get_OneSNP_Geno scans) still read packed_byte
	// and would crash. Set SAIGE_GPU_RELEASE_HOST=1 for GPU-only runs
	// that need the host-side memory for other allocations. Note: this
	// does NOT reduce PEAK RSS, because the peak is set during BED load +
	// GPU upload (both buffers live simultaneously).
	// Only legal for tiers 3/4: those upload once and never dereference the
	// host buffer again. Tiers 1/2 re-read packed_flat_ on every matvec, so
	// releasing it there would be a use-after-free.
	{
		static int  s_success_count = 0;
		static bool s_released      = false;
		static const bool s_release_enabled = [](){
			const char* v = std::getenv("SAIGE_GPU_RELEASE_HOST");
			return v && std::string(v) == "1";
		}();
		// SAIGE_GPU_VERIFY=1: run the CPU path too and print the rel-L2 gap.
		// Mirrors SAIGE_GEMV_VERIFY / SAIGE_AVX2_VERIFY. Roughly halves
		// throughput, so it is a debugging switch, not a production one.
		static const bool s_gpu_verify = [](){
			const char* v = std::getenv("SAIGE_GPU_VERIFY");
			return v && std::string(v) == "1";
		}();
		if (s_gpu_handle) {
			arma::fvec out((arma::uword)geno.getNnomissing());
			if (saige::gpu::matvec(s_gpu_handle, bVec.memptr(), out.memptr())) {
				if (s_gpu_verify) {
					g_gpu_in_verify = true;
					arma::fvec ref = parallelCrossProd(bVec);   // CPU reference
					g_gpu_in_verify = false;
					const double rn = arma::norm(ref);
					const double rel = arma::norm(out - ref) / std::max(1e-30, rn);
					static long s_vcalls = 0;
					std::cout << "[GPU VERIFY call#" << (++s_vcalls)
					          << "] |ref|=" << rn
					          << "  max_abs=" << arma::abs(out - ref).max()
					          << "  rel_l2=" << rel
					          << (rel < 1e-5 ? "  OK" : "  **ABOVE 1e-5**")
					          << std::endl;
				}
				if (s_release_enabled && !s_released && ++s_success_count >= 4 &&
				    saige::gpu::tier(s_gpu_handle) >= 3 && geno.use_packed_flat_) {
					geno.release_packed_flat_host();
					s_released = true;
				}
				return out;
			}
			// one-time failure → fall through to CPU on this and all future calls
			gpu_mark_failed("parallelCrossProd");
		}
	}

	int Msub_mafge1perc = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();

	// Phase-1: optional psi-v wall-time telemetry (SAIGE_TIME_PSIV=1).
	// Prints cumulative call count / total / mean every 50 calls.
	static const bool s_time_psiv = [](){
		const char* v = std::getenv("SAIGE_TIME_PSIV");
		return v && std::string(v) == "1"; }();
	struct PsivTimer {
		bool on; double t0 = 0;
		static double& total() { static double t = 0; return t; }
		static long&   calls() { static long c = 0; return c; }
		PsivTimer(bool enabled) : on(enabled) { if (on) t0 = get_wall_time(); }
		~PsivTimer() {
			if (!on) return;
			total() += get_wall_time() - t0;
			if ((++calls()) % 50 == 0)
				std::cout << "[PSIV] calls=" << calls() << " total=" << total()
				          << "s avg=" << 1e3 * total() / calls() << "ms"
				          << std::endl;
		}
	} psiv_timer(s_time_psiv);

	// DEBUG: Detailed diagnostics on first call
	static int debug_call_count = 0;
	if(debug_call_count == 0) {
		std::cout << "\n=== DEBUG parallelCrossProd (first call) ===" << std::endl;
		std::cout << "Msub_mafge1perc (expected markers): " << Msub_mafge1perc << std::endl;
		std::cout << "Nnomissing: " << geno.getNnomissing() << std::endl;
		std::cout << "N_fam: " << geno.getN() << std::endl;
		std::cout << "bVec size: " << bVec.n_elem << std::endl;
		std::cout << "bVec[0:5]: " << bVec[0] << " " << bVec[1] << " " << bVec[2]
		          << " " << bVec[3] << " " << bVec[4] << std::endl;

		arma::fvec vec_debug;
		// Check first 3 markers: freq, invstd, genotype distribution, stdgeno values
		for(int m = 0; m < 3; m++) {
			geno.Get_OneSNP_StdGeno(m, &vec_debug);
			float freq_m = geno.alleleFreqVec[m];
			float invstd_m = geno.invstdvVec[m];
			// Count genotype distribution from stdgeno values
			float std0 = (0 - 2*freq_m) * invstd_m;  // geno=0
			float std1 = (1 - 2*freq_m) * invstd_m;  // geno=1
			float std2 = (2 - 2*freq_m) * invstd_m;  // geno=2
			int cnt0=0, cnt1=0, cnt2=0, cntother=0;
			for(size_t s=0; s < vec_debug.n_elem; s++) {
				float v = vec_debug[s];
				if(std::fabs(v - std0) < 1e-4) cnt0++;
				else if(std::fabs(v - std1) < 1e-4) cnt1++;
				else if(std::fabs(v - std2) < 1e-4) cnt2++;
				else cntother++;
			}
			float dot_val = arma::dot(vec_debug, bVec);
			std::cout << "Marker " << m << ": freq=" << freq_m << " invstd=" << invstd_m
			          << " geno_counts=[" << cnt0 << "," << cnt1 << "," << cnt2 << "] other=" << cntother
			          << " |vec|=" << arma::norm(vec_debug)
			          << " dot=" << dot_val << std::endl;
			std::cout << "  StdGeno[0:5]: " << vec_debug[0] << " " << vec_debug[1] << " "
			          << vec_debug[2] << " " << vec_debug[3] << " " << vec_debug[4] << std::endl;
			std::cout << "  expected_std: g0=" << std0 << " g1=" << std1 << " g2=" << std2 << std::endl;
		}
		// Also check marker stats: mean and variance of stdgeno should be ~0 and ~1
		geno.Get_OneSNP_StdGeno(0, &vec_debug);
		std::cout << "Marker 0: mean(stdgeno)=" << arma::mean(vec_debug)
		          << " var(stdgeno)=" << arma::var(vec_debug) << std::endl;
		debug_call_count++;
	}

	CorssProd CorssProd(bVec);

	// ------------------------------------------------------------------
	// Phase-1 AVX2 fused decode (AVX2_KERNEL_PLAN.md). Default ON when the
	// binary was built with AVX2+FMA and per-marker rows are contiguous.
	// Env toggles:
	//   SAIGE_NO_AVX2=1         -> original scalar path
	//   SAIGE_AVX2_PASS1_ONLY=1 -> AVX2 val1 only, scalar axpy (validation)
	//   SAIGE_AVX2_VERIFY=1     -> also run scalar path, print rel-L2 diff
	// ------------------------------------------------------------------
	int avx2_mode = 0;
	arma::fvec xperm;
#if SAIGE_AVX2_KERNEL_AVAILABLE
	static const bool s_no_avx2 = [](){
		const char* v = std::getenv("SAIGE_NO_AVX2");
		return v && std::string(v) == "1"; }();
	static const bool s_pass1_only = [](){
		const char* v = std::getenv("SAIGE_AVX2_PASS1_ONLY");
		return v && std::string(v) == "1"; }();
	static const bool s_avx2_verify = [](){
		const char* v = std::getenv("SAIGE_AVX2_VERIFY");
		return v && std::string(v) == "1"; }();
	if (!s_no_avx2 && geno.packed_rows_contiguous()) {
		const arma::uword N = geno.getNnomissing();
		xperm.set_size(N);
		saige_avx2::permute_fwd(bVec.memptr(), xperm.memptr(), (std::size_t)N);
		double Sx = 0.0;
		const float* bp = bVec.memptr();
		for (arma::uword i = 0; i < N; ++i) Sx += bp[i];
		avx2_mode = s_pass1_only ? 1 : 2;
		CorssProd.m_avx2_mode = avx2_mode;
		CorssProd.m_xperm     = xperm.memptr();
		CorssProd.m_Sx        = (float)Sx;
		static bool s_announced = false;
		if (!s_announced) {
			std::cout << "[parallelCrossProd] AVX2 fused-decode kernel enabled (mode "
			          << avx2_mode << ")." << std::endl;
			s_announced = true;
		}
	}
#endif

  	parallelReduce(0, Msub_mafge1perc, CorssProd);

#if SAIGE_AVX2_KERNEL_AVAILABLE
	if (avx2_mode == 2) {
		const arma::uword N = geno.getNnomissing();
		arma::fvec out(N);
		saige_avx2::unpermute_sub(CorssProd.m_bout.memptr(), out.memptr(),
		                          (std::size_t)N, CorssProd.m_Coffset);
		out /= static_cast<float>(CorssProd.Msub_mafge1perc);
		if (s_avx2_verify) {
			struct CorssProd ref(bVec);
			parallelReduce(0, Msub_mafge1perc, ref);
			arma::fvec ref_out =
			    ref.m_bout / static_cast<float>(ref.Msub_mafge1perc);
			const double diff = arma::norm(out - ref_out) /
			                    std::max(1e-30, (double)arma::norm(ref_out));
			std::cout << "[AVX2 VERIFY] psi-v rel-L2 diff vs scalar = " << diff
			          << (diff < 1e-5 ? "  OK" : "  **ABOVE 1e-5**") << std::endl;
		}
		return out;
	}
	if (avx2_mode == 1 && s_avx2_verify) {
		struct CorssProd ref(bVec);
		parallelReduce(0, Msub_mafge1perc, ref);
		const double diff = arma::norm(CorssProd.m_bout - ref.m_bout) /
		                    std::max(1e-30, (double)arma::norm(ref.m_bout));
		std::cout << "[AVX2 VERIFY] pass1-only rel-L2 diff vs scalar = " << diff
		          << (diff < 1e-5 ? "  OK" : "  **ABOVE 1e-5**") << std::endl;
	}
#endif

	if(debug_call_count == 1) {
		// Dump first 3 markers' full genotype vectors to file for R comparison
		{
			std::string dump_path = "/data/ukb_cpp_test/output/cpp_stdgeno_dump.csv";
			std::ofstream dump(dump_path);
			if (dump.is_open()) {
				arma::fvec vec_dump;
				dump << "sample_idx";
				for (int m = 0; m < 3; m++) dump << ",marker" << m << "_stdgeno,marker" << m << "_rawgeno";
				dump << "\n";
				// Get raw genotypes via Get_OneSNP_Geno
				arma::fvec svec0, svec1, svec2;
				geno.Get_OneSNP_StdGeno(0, &svec0);
				geno.Get_OneSNP_StdGeno(1, &svec1);
				geno.Get_OneSNP_StdGeno(2, &svec2);
				arma::ivec* raw0 = geno.Get_OneSNP_Geno(0);
				arma::ivec* raw1 = geno.Get_OneSNP_Geno(1);
				arma::ivec* raw2 = geno.Get_OneSNP_Geno(2);
				int n_dump = std::min((int)svec0.n_elem, 100);  // first 100 samples
				for (int s = 0; s < n_dump; s++) {
					dump << s
					     << "," << svec0[s] << "," << (*raw0)[s]
					     << "," << svec1[s] << "," << (*raw1)[s]
					     << "," << svec2[s] << "," << (*raw2)[s]
					     << "\n";
				}
				dump.close();
				std::cout << "Dumped first 3 markers x 100 samples to: " << dump_path << std::endl;
			}
		}
		std::cout << "CorssProd.Msub_mafge1perc (actual markers): " << CorssProd.Msub_mafge1perc << std::endl;
		std::cout << "bout[0:5] (before /M, float): " << CorssProd.m_bout[0] << " " << CorssProd.m_bout[1]
		          << " " << CorssProd.m_bout[2] << " " << CorssProd.m_bout[3] << " " << CorssProd.m_bout[4] << std::endl;
		arma::fvec result = CorssProd.m_bout / static_cast<float>(CorssProd.Msub_mafge1perc);
		std::cout << "result[0:5] (Ku = bout/M): " << result[0] << " " << result[1]
		          << " " << result[2] << " " << result[3] << " " << result[4] << std::endl;
		float corr = arma::as_scalar(arma::cor(result, bVec));
		std::cout << "cor(Ku, u) = " << corr << " (should be <1 if K has off-diagonal structure)" << std::endl;
		std::cout << "=== END DEBUG parallelCrossProd ===" << std::endl << std::endl;
		debug_call_count++;
	}

  	return CorssProd.m_bout / static_cast<float>(CorssProd.Msub_mafge1perc);
}


// ----------------------------------------------------------------------------
// Blocked-GEMV K·u prototype.
//
// Replaces the SNP-by-SNP `dot + axpy` loop in `parallelCrossProd` with two
// BLAS GEMVs per block of B markers. Layout: a single reusable arma::fmat
// `Z_block` of shape (N, B) where each column is a standardized marker vector
// (sample-contiguous in arma's column-major storage) — that gives:
//   - contiguous writes when decoding a marker into a column
//   - efficient GEMV in both directions
//
// Math: for a chunk of markers m ∈ [m0, m0+b):
//   tmp = Z_block(:, 0:b).t() * bVec   // length b
//   out += Z_block(:, 0:b)     * tmp   // length N
//
// Final result is `out / Msub_mafge1perc`, matching the original
// `parallelCrossProd` normalization. Threading: the decode loop is serial,
// BLAS handles intra-GEMV parallelism. Set OPENBLAS_NUM_THREADS / MKL etc.
// via the existing `nthreads` plumbing; this function does not touch BLAS
// thread state.
//
// Toggles (env, set by main.cpp from cfg / CLI):
//   SAIGE_USE_BLOCKED_GEMV=1   -- enable
//   SAIGE_GEMV_BLOCK_SIZE=B    -- block size (default 128)
//   SAIGE_GEMV_VERIFY=1        -- A/B compare against original on every call
//
arma::fvec parallelCrossProd_blocked(arma::fcolvec & bVec) {
	const int M = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
	const int N = (int)geno.getNnomissing();

	// Block size from env (default 128). Clamp into [16, 4096].
	int B = 128;
	if (const char* bs = std::getenv("SAIGE_GEMV_BLOCK_SIZE")) {
		int v = std::atoi(bs);
		if (v >= 16 && v <= 4096) B = v;
	}
	if (B > M) B = M;

	// Reusable buffers across calls. Layout: Z_block(N, B) — columns are markers.
	static arma::fmat Z_block;
	static arma::fvec tmp_buf;
	if ((int)Z_block.n_rows != N || (int)Z_block.n_cols != B) {
		Z_block.set_size(N, B);
	}
	if ((int)tmp_buf.n_elem != B) {
		tmp_buf.set_size(B);
	}

	// Per-marker decode scratch. Reused across markers — Get_OneSNP_StdGeno
	// will resize this to N on first use, then in-place writes.
	arma::fvec vec_scratch;

	arma::fvec out(N, arma::fill::zeros);

	for (int m0 = 0; m0 < M; m0 += B) {
		const int b_eff = std::min(B, M - m0);

		// 1) Decode b_eff markers into columns of Z_block.
		// Each column is sample-contiguous (column-major) so writes are linear.
		for (int k = 0; k < b_eff; ++k) {
			geno.Get_OneSNP_StdGeno((size_t)(m0 + k), &vec_scratch);
			// Copy decoded marker into column k of Z_block.
			std::memcpy(Z_block.colptr(k),
			            vec_scratch.memptr(),
			            sizeof(float) * (size_t)N);
		}

		// 2) Two GEMVs. Use submatrix view when the last block is partial.
		if (b_eff == B) {
			arma::fvec tmp = Z_block.t() * bVec;       // (B x N) * (N) -> (B)
			out += Z_block * tmp;                      // (N x B) * (B) -> (N)
		} else {
			arma::fmat Zsub = Z_block.cols(0, b_eff - 1);
			arma::fvec tmp  = Zsub.t() * bVec;
			out += Zsub * tmp;
		}
	}

	return out / static_cast<float>(M);
}


// ----------------------------------------------------------------------------
// Phase-2 multi-RHS ψ·B (block-PCG support). Computes ψ·B for a whole N×k
// matrix of right-hand sides while streaming the packed genotype matrix
// (the DRAM-bound resource) only twice — the same traffic as ONE
// single-column ψv call.
//
// Structure per TBB worker (marker range [begin,end)):
//   sweep 1 (dot):  raw[m][j] = Σ_i g[m,i]·Xp[i,j], sample-blocked so the
//                   Xp slice (k×SBS floats) stays L2-resident across the
//                   whole marker range (unblocked, the N×k RHS matrix is
//                   re-read from L3/DRAM per marker — measured to erase the
//                   entire batching win at k≳8).
//   sval[m][j]     = invstd²·(raw − 2f·Sx[j]) per rank-one identity
//   sweep 2 (axpy): bout[:,j] += sval[m][j]·g[m,:], sample-blocked the same
//                   way so the bout slice stays L2-resident.
// bout accumulates in PERM order without the uniform −2f·sval offset
// (tracked in m_Coffset per column); the caller un-permutes and subtracts.
//
// Scalar fallback (mode 0, SAIGE_NO_AVX2=1 or non-contiguous storage):
// decode each marker once via Get_OneSNP_StdGeno, then k dot/axpy pairs.
// ----------------------------------------------------------------------------
struct CorssProdMat : public Worker
{
	const arma::fmat& m_Bmat;   // N×k RHS, natural order (scalar mode)
	unsigned int m_N;
	unsigned int m_k;

	arma::fmat  m_bout;         // N×k accumulator (perm order in avx2 mode)
	arma::fvec  m_Coffset;      // k per-column uniform offsets (avx2 mode)
	int Msub_mafge1perc;

	int          m_avx2_mode = 0;        // 0 scalar, 2 avx2 two-sweep
	const float* m_xperm     = nullptr;  // N×k perm-order RHS (ldx = N)
	const float* m_Sx        = nullptr;  // k column sums of Bmat
	std::size_t  m_sbs       = 4096;     // sample-block size (multiple of 128)

	CorssProdMat(const arma::fmat& B)
		: m_Bmat(B)
	{
		m_N = geno.getNnomissing();
		m_k = B.n_cols;
		m_bout.zeros(m_N, m_k);
		m_Coffset.zeros(m_k);
		Msub_mafge1perc = 0;
	}
	CorssProdMat(const CorssProdMat& other, Split)
		: m_Bmat(other.m_Bmat)
	{
		m_N = other.m_N;
		m_k = other.m_k;
		m_bout.zeros(m_N, m_k);
		m_Coffset.zeros(m_k);
		Msub_mafge1perc = 0;
		m_avx2_mode = other.m_avx2_mode;
		m_xperm     = other.m_xperm;
		m_Sx        = other.m_Sx;
		m_sbs       = other.m_sbs;
	}

	void operator()(std::size_t begin, std::size_t end) {
#if SAIGE_AVX2_KERNEL_AVAILABLE
		if (m_avx2_mode != 0) {
			const std::size_t nm = end - begin;
			const std::size_t k  = m_k;
			const std::size_t N  = m_N;
			// raw dot accumulators + per-marker scale factors, row-major [m][j]
			std::vector<float> raws(nm * k, 0.0f);
			std::vector<float> tmp(k);

			// sweep 1: sample-blocked dot products
			for (std::size_t sb = 0; sb < N; sb += m_sbs) {
				const std::size_t len = std::min(m_sbs, N - sb);
				const float* xp_sb = m_xperm + sb;  // col j at + j*N
				for (std::size_t m = begin; m < end; ++m) {
					const unsigned char* row =
					    geno.packed_row_ptr(m) + sb / 4;
					saige_avx2::pass1_sum_gx_multi(row, xp_sb, N, k, len,
					                               tmp.data());
					float* r = raws.data() + (m - begin) * k;
					for (std::size_t j = 0; j < k; ++j) r[j] += tmp[j];
				}
			}

			// rank-one identity: raw -> sval; accumulate Coffset
			for (std::size_t m = begin; m < end; ++m) {
				const float f = geno.alleleFreqVec[m];
				const float d = geno.invstdvVec[m];
				float* r = raws.data() + (m - begin) * k;
				for (std::size_t j = 0; j < k; ++j) {
					const float val1 = d * (r[j] - 2.0f * f * m_Sx[j]);
					const float sval = val1 * d;
					r[j] = sval;                        // reuse as sval[m][j]
					m_Coffset[j] += 2.0f * f * sval;
				}
			}

			// sweep 2: sample-blocked axpy into perm-order bout
			float* boutp = m_bout.memptr();             // col j at + j*N
			for (std::size_t sb = 0; sb < N; sb += m_sbs) {
				const std::size_t len = std::min(m_sbs, N - sb);
				for (std::size_t m = begin; m < end; ++m) {
					const unsigned char* row =
					    geno.packed_row_ptr(m) + sb / 4;
					saige_avx2::pass2_axpy_multi(
					    row, raws.data() + (m - begin) * k, boutp + sb, N, k,
					    len);
				}
			}
			Msub_mafge1perc += (int)nm;
			return;
		}
#endif
		// scalar fallback: decode once per marker, k dot/axpy pairs
		arma::fvec vec;
		for (std::size_t m = begin; m < end; ++m) {
			geno.Get_OneSNP_StdGeno(m, &vec);
			for (unsigned int j = 0; j < m_k; ++j) {
				const float val1 = arma::dot(vec, m_Bmat.col(j));
				m_bout.col(j) += val1 * vec;
			}
			Msub_mafge1perc += 1;
		}
	}

	void join(const CorssProdMat& rhs) {
		m_bout += rhs.m_bout;
		m_Coffset += rhs.m_Coffset;
		Msub_mafge1perc += rhs.Msub_mafge1perc;
	}
};

// ψ·B for an N×k RHS matrix (dense-GRM path). Multi-column analogue of
// parallelCrossProd; same normalization (division by Msub_mafge1perc).
arma::fmat parallelCrossProdMat(const arma::fmat& Bmat) {
	const int Msub_mafge1perc =
	    geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
	const arma::uword N = geno.getNnomissing();
	const arma::uword k = Bmat.n_cols;

	// optional ψ·B wall-time telemetry (SAIGE_TIME_PSIV=1), separate
	// accumulator from the single-column [PSIV] one.
	static const bool s_time_psiv = [](){
		const char* v = std::getenv("SAIGE_TIME_PSIV");
		return v && std::string(v) == "1"; }();
	struct PsivMatTimer {
		bool on; double t0 = 0; arma::uword k;
		static double& total() { static double t = 0; return t; }
		static long&   calls() { static long c = 0; return c; }
		static long&   cols()  { static long c = 0; return c; }
		PsivMatTimer(bool enabled, arma::uword kk) : on(enabled), k(kk) {
			if (on) t0 = get_wall_time(); }
		~PsivMatTimer() {
			if (!on) return;
			total() += get_wall_time() - t0;
			cols()  += (long)k;
			if ((++calls()) % 10 == 0)
				std::cout << "[PSIVMAT] calls=" << calls() << " cols=" << cols()
				          << " total=" << total() << "s avg/col="
				          << 1e3 * total() / cols() << "ms" << std::endl;
		}
	} psivmat_timer(s_time_psiv, k);

	// G6: batch GPU dispatch, sharing the single-column handle. Only tier 4
	// has a multi-RHS kernel; every other tier answers "not available" and we
	// stay on the CPU rather than looping matvec() k times (that would re-read
	// the whole packed matrix once per column, which is the very cost the
	// batch kernel exists to avoid).
	if (saige::gpu::Handle* h = gpu_handle_or_null()) {
		if (saige::gpu::matvec_mat_available(h)) {
			arma::fmat out(N, k);
			if (saige::gpu::matvec_mat(h, Bmat.memptr(), (int)k, out.memptr()))
				return out;
			gpu_mark_failed("parallelCrossProdMat");
		}
	}

	CorssProdMat worker(Bmat);

	arma::fmat xperm;
	arma::fvec Sx;
#if SAIGE_AVX2_KERNEL_AVAILABLE
	static const bool s_no_avx2 = [](){
		const char* v = std::getenv("SAIGE_NO_AVX2");
		return v && std::string(v) == "1"; }();
	if (!s_no_avx2 && geno.packed_rows_contiguous()) {
		xperm.set_size(N, k);
		Sx.set_size(k);
		for (arma::uword j = 0; j < k; ++j) {
			saige_avx2::permute_fwd(Bmat.colptr(j), xperm.colptr(j),
			                        (std::size_t)N);
			double s = 0.0;
			const float* bp = Bmat.colptr(j);
			for (arma::uword i = 0; i < N; ++i) s += bp[i];
			Sx(j) = (float)s;
		}
		worker.m_avx2_mode = 2;
		worker.m_xperm     = xperm.memptr();
		worker.m_Sx        = Sx.memptr();
		// sample-block size: keep k×SBS floats (the RHS slice in sweep 1,
		// the bout slice in sweep 2) within ~96 KB of the 256 KB L2.
		// SAIGE_BLOCKPCG_SBS overrides (0 = no blocking / whole row).
		static const long s_sbs_env = [](){
			const char* v = std::getenv("SAIGE_BLOCKPCG_SBS");
			return v ? std::atol(v) : -1L; }();
		std::size_t sbs;
		if (s_sbs_env == 0) {
			sbs = (std::size_t)N;  // unblocked
		} else if (s_sbs_env > 0) {
			sbs = ((std::size_t)s_sbs_env) & ~(std::size_t)127;
		} else {
			sbs = (std::size_t)(24576 / std::max<arma::uword>(k, 1));
			sbs &= ~(std::size_t)127;                   // multiple of 128
		}
		worker.m_sbs = std::max<std::size_t>(sbs, 512); // >= 4 superblocks
		static bool s_announced = false;
		if (!s_announced) {
			std::cout << "[parallelCrossProdMat] AVX2 multi-RHS kernel enabled"
			          << " (sbs=" << worker.m_sbs << ")." << std::endl;
			s_announced = true;
		}
	}
#endif

	// coarse grain: sample-blocking amortizes the RHS slice over the marker
	// range, so avoid tiny TBB chunks.
	const std::size_t grain =
	    std::max<std::size_t>(64, (std::size_t)Msub_mafge1perc / 32);
	parallelReduce(0, Msub_mafge1perc, worker, grain);

#if SAIGE_AVX2_KERNEL_AVAILABLE
	if (worker.m_avx2_mode == 2) {
		arma::fmat out(N, k);
		for (arma::uword j = 0; j < k; ++j)
			saige_avx2::unpermute_sub(worker.m_bout.colptr(j), out.colptr(j),
			                          (std::size_t)N, worker.m_Coffset(j));
		out /= static_cast<float>(worker.Msub_mafge1perc);
		return out;
	}
#endif
	return worker.m_bout / static_cast<float>(worker.Msub_mafge1perc);
}


// REMOVED: innerProductFun() - use getInnerProd() from src/UTIL.cpp instead


/*
// INTERNAL: LOCO-specific parallel cross product computation
arma::fvec parallelCrossProd_LOCO_2(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int Msub = geno.getMsub();
	int M = geno.getM();
	
        CorssProd_LOCO CorssProd_LOCO(bVec);

  // call paralleReduce to start the work
        parallelReduce(0, M, CorssProd_LOCO);

  // return the computed product
	//cout << "Msub: " << Msub << endl;

	for(int i=0; i<10; ++i)
        {
		std::cout << (CorssProd_LOCO.m_bout)[i] << ' ';
        }
	std::cout << std::endl;
	std::cout << "CorssProd_LOCO.m_Msub_mafge1perc: " << CorssProd_LOCO.m_Msub_mafge1perc << std::endl;
        //return CorssProd_LOCO.m_bout/Msub;
	return CorssProd_LOCO.m_bout/(CorssProd_LOCO.m_Msub_mafge1perc);
}
*/



// INTERNAL: Full parallel cross product computation with marker count
arma::fvec parallelCrossProd_full(arma::fcolvec & bVec, int & markerNum) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int M = geno.getM();
	//
        int Msub_mafge1perc = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM(); 
        CorssProd CorssProd(bVec);

        //std::cout << "Msub_mafge1perc ok  " << Msub_mafge1perc << std::endl;        
  // call paralleReduce to start the work
        parallelReduce(0, Msub_mafge1perc, CorssProd);
        markerNum = CorssProd.Msub_mafge1perc;
        //std::cout << "markerNum " << markerNum << std::endl;        

        //cout << "print test; M: " << M << endl;
        //for(int i=0; i<10; ++i)
        //{
        //        cout << (CorssProd.m_bout)[i] << ' ' << endl;
        //        cout << bVec[i] << ' ' << endl;
        //        cout << (CorssProd.m_bout/M)[i] << ' ' << endl;
        //}
        ////cout << endl;
  // return the computed product
        //std::cout << "number of markers with maf ge " << minMAFtoConstructGRM << " is " << CorssProd.Msub_mafge1perc << std::endl;
        return arma::conv_to<arma::fvec>::from(CorssProd.m_bout);
}


// INTERNAL: LOCO parallel cross product computation
arma::fvec parallelCrossProd_LOCO(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int Msub = geno.getMsub();
	//int M = geno.getM();
        int numberMarker_full = 0;
        arma::fvec outvec = parallelCrossProd_full(bVec, numberMarker_full);

        //CorssProd_LOCO CorssProd_LOCO(bVec);
	CorssProd CorssProd(bVec);
  // call paralleReduce to start the work
	int startIndex = geno.getStartIndex();
        int endIndex = geno.getEndIndex();

	parallelReduce(startIndex, endIndex+1, CorssProd);



	outvec = outvec - arma::conv_to<arma::fvec>::from(CorssProd.m_bout);

	/*
	for(int i=0; i<10; ++i)
        {
                std::cout << (outvec)[i] << ' ';
        }
        std::cout << std::endl;
*/

	int markerNum = numberMarker_full - CorssProd.Msub_mafge1perc;
	//std::cout << "markerNum: " << markerNum << std::endl; 
      	// return the computed product
	//cout << "Msub: " << Msub << endl;
        //for(int i=0; i<100; ++i)
        //{
        //	cout << (CorssProd_LOCO.m_bout/Msub)[i] << ' ';
        //}
	//cout << endl;
        //return CorssProd_LOCO.m_bout/Msub;
        //return CorssProd_LOCO.m_bout/(CorssProd_LOCO.m_Msub_mafge1perc);
	return outvec/markerNum;
}


arma::umat locationMat;
arma::vec valueVec;
int dimNum = 0;

// INTERNAL: Utility function for setting up sparse genetic relationship matrix
void setupSparseGRM(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    // sparse x sparse -> sparse
    //arma::sp_mat result(a);
    //int r = a.n_rows;
        locationMat.zeros(2,r);
        valueVec.zeros(r);

    locationMat = locationMatinR;
    valueVec = valueVecinR;
    dimNum = r;

    std::cout << locationMat.n_rows << " locationMat.n_rows " << std::endl;
    std::cout << locationMat.n_cols << " locationMat.n_cols " << std::endl;
    std::cout << valueVec.n_elem << " valueVec.n_elem " << std::endl;
    
    //for(size_t i=0; i< 10; i++){
    //    std::cout << valueVec(i) << std::endl;
    //    std::cout << locationMat(0,i) << std::endl;
    //    std::cout << locationMat(1,i) << std::endl;
    //}

    //arma::vec y = arma::linspace<arma::vec>(0, 5, r);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    //arma::vec x = arma::spsolve( result, y );

    //return x;
}

// Build sparse GRM (uses geno + internal kernels)
void build_sparse_grm_in_place(double relatedness_cutoff,
                               double min_maf, double max_miss)
{
  std::cout << "[DEBUG build_sparse_grm_in_place] START" << std::endl;
  std::cout << "[DEBUG] relatedness_cutoff=" << relatedness_cutoff
            << " min_maf=" << min_maf << " max_miss=" << max_miss << std::endl;

  setminMAFforGRM(static_cast<float>(min_maf));
  setmaxMissingRateforGRM(static_cast<float>(max_miss));
  setRelatednessCutoff(static_cast<float>(relatedness_cutoff));

  // choose markers (MAF QC)
  std::cout << "[DEBUG] Step 1: getQCdMarkerIndex()..." << std::endl;
  std::vector<bool> keep = getQCdMarkerIndex();
  int nKeep = std::count(keep.begin(), keep.end(), true);
  std::cout << "[DEBUG] Markers passing QC: " << nKeep << " / " << keep.size() << std::endl;

  // FIX BUG #3: Use COMPACTED indices [0, 1, 2, ..., nKeep-1]
  // NOT original indices, because genoVecofPointers/alleleFreqVec/invstdvVec
  // are stored at compacted indices during loading (SNPIdx_new)
  arma::ivec sub; sub.set_size(nKeep);
  for (int i = 0; i < nKeep; ++i) sub[i] = i;  // compacted indices

  std::cout << "[DEBUG] Step 2: setSubMarkerIndex()..." << std::endl;
  setSubMarkerIndex(sub);
  std::cout << "[DEBUG] subMarkerIndex.n_elem=" << sub.n_elem
            << " Nnomissing=" << geno.getNnomissing() << std::endl;
  std::cout << "[DEBUG] stdGenoMultiMarkersMat dimensions after setSubMarkerIndex: "
            << geno.stdGenoMultiMarkersMat.n_rows << " x " << geno.stdGenoMultiMarkersMat.n_cols << std::endl;

  // std-geno blocks for pair screening
  std::cout << "[DEBUG] Step 3: Get_MultiMarkersBySample_StdGeno_Mat()..." << std::endl;
  std::cout << "[DEBUG] genoVecofPointers.size()=" << geno.genoVecofPointers.size()
            << " numMarkersofEachArray=" << geno.numMarkersofEachArray
            << " m_size_of_esi=" << geno.m_size_of_esi << std::endl;
  Get_MultiMarkersBySample_StdGeno_Mat();
  std::cout << "[DEBUG] After Get_MultiMarkersBySample_StdGeno_Mat: matrix "
            << geno.stdGenoMultiMarkersMat.n_rows << " x " << geno.stdGenoMultiMarkersMat.n_cols << std::endl;

  // related pairs (FIX: kinship values are now computed and stored in findIndiceRelatedSample)
  std::cout << "[DEBUG] Step 4: findIndiceRelatedSample()..." << std::endl;
  findIndiceRelatedSample();
  std::cout << "[DEBUG] findIndiceRelatedSample() COMPLETED" << std::endl;
  int npairs = (int)geno.indiceVec.size();
  std::cout << "[DEBUG] npairs (related sample pairs found)=" << npairs << std::endl;
  std::cout << "[DEBUG] kinValueVecSparse.size()=" << geno.kinValueVecSparse.size() << std::endl;

  // FIX: Use kinship values already computed in findIndiceRelatedSample
  // (parallelcalsparseGRM was using m_OneSNP_Geno which isn't set up here)

  // assemble COO + add diagonal = 1.0 (store both triangles for symmetry)
  const int n  = geno.getNnomissing();
  const int nz = 2*npairs + n;  // both (i,j) and (j,i) + diagonal
  arma::umat loc(2, nz); arma::vec val(nz);
  for (int t=0; t<npairs; ++t) {
    loc(0,t) = (unsigned)geno.indiceVec[t].first;
    loc(1,t) = (unsigned)geno.indiceVec[t].second;
    val(t)   = (double)geno.kinValueVecSparse[t];
    // Store symmetric entry (j,i)
    loc(0, npairs+t) = (unsigned)geno.indiceVec[t].second;
    loc(1, npairs+t) = (unsigned)geno.indiceVec[t].first;
    val(npairs+t)    = (double)geno.kinValueVecSparse[t];
  }
  for (int i=0; i<n; ++i) {
    loc(0, 2*npairs+i) = i;
    loc(1, 2*npairs+i) = i;
    val(   2*npairs+i) = 1.0;
  }
  setupSparseGRM(n, loc, val);
}

// Exporters
arma::umat export_sparse_grm_locations() { return locationMat; }
arma::vec  export_sparse_grm_values()    { return valueVec; }
int        export_sparse_grm_dim()       { return dimNum; }
int        export_sparse_grm_nnz()       { return (int)valueVec.n_elem; }



bool isUsePrecondM = false;
bool isUseSparseSigmaforInitTau = false;
bool isUseSparseSigmaforModelFitting = false;
bool isUsePCGwithSparseSigma = false;

// Status getter (returns the internal flag)
bool get_isUseSparseSigmaforModelFitting() { return isUseSparseSigmaforModelFitting; }

// INTERNAL: Compute cross product with kinship matrix  
arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec){
       arma::fvec crossProdVec;

    // Debug: print input/output for first few calls (compile with -DSAIGE_DEBUG_IO to enable)
    static int crossprod_call_count = 0;
#ifdef SAIGE_DEBUG_IO
    if (crossprod_call_count < 3) {
        std::cout << "[getCrossprodMatAndKin call #" << crossprod_call_count << "]"
                  << " input |u|=" << arma::norm(bVec)
                  << " u[0:5]=" << bVec[0] << " " << bVec[1] << " " << bVec[2] << " " << bVec[3] << " " << bVec[4]
                  << " n=" << bVec.n_elem
                  << " sparse=" << (isUseSparseSigmaforInitTau | isUseSparseSigmaforModelFitting)
                  << std::endl;
    }
#endif

    if(isUseSparseSigmaforInitTau | isUseSparseSigmaforModelFitting){
        //cout << "use sparse kinship to estimate initial tau and for getCrossprodMatAndKin" <<  endl;


	arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
	arma::vec x = result * arma::conv_to<arma::dcolvec>::from(bVec);


//double wall3in = get_wall_time();
// double cpu3in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;


    crossProdVec = arma::conv_to<arma::fvec>::from(x);

   }else{
        // Dispatch: blocked GEMV path (CPU only, prototype) vs. original SNP-loop.
        // Sparse-Sigma path above takes precedence so the flag is only consulted
        // here, on the dense GRM branch where parallelCrossProd is the inner loop.
        static const bool s_use_blocked = [](){
            const char* v = std::getenv("SAIGE_USE_BLOCKED_GEMV");
            return v && std::string(v) == "1";
        }();
        static const bool s_verify = [](){
            const char* v = std::getenv("SAIGE_GEMV_VERIFY");
            return v && std::string(v) == "1";
        }();

        if (s_verify) {
            // Run both paths and report rel/abs error. Original path is the reference.
            arma::fvec ref = parallelCrossProd(bVec);
            arma::fvec alt = parallelCrossProd_blocked(bVec);
            arma::fvec diff = alt - ref;
            float ref_norm = arma::norm(ref);
            float abs_max  = arma::abs(diff).max();
            float rel_l2   = (ref_norm > 0.0f) ? (arma::norm(diff) / ref_norm) : 0.0f;
            std::cout << "[GEMV_VERIFY call#" << crossprod_call_count
                      << "] |ref|=" << ref_norm
                      << "  max_abs_err=" << abs_max
                      << "  rel_l2_err=" << rel_l2
                      << std::endl;
            // Use the blocked result when the flag is on so the rest of the
            // solver consumes the new path; otherwise use the reference.
            crossProdVec = s_use_blocked ? alt : ref;
        } else if (s_use_blocked) {
            crossProdVec = parallelCrossProd_blocked(bVec);
        } else {
            crossProdVec = parallelCrossProd(bVec);
        }
   }

#ifdef SAIGE_DEBUG_IO
    if (crossprod_call_count < 3) {
        std::cout << "[C++ getCrossprodMatAndKin call #" << crossprod_call_count << "]"
                  << " |u|=" << arma::norm(bVec)
                  << " |Au|=" << arma::norm(crossProdVec)
                  << " u[0:5]=" << bVec[0] << " " << bVec[1] << " " << bVec[2] << " " << bVec[3] << " " << bVec[4]
                  << " Au[0:5]=" << crossProdVec[0] << " " << crossProdVec[1] << " " << crossProdVec[2] << " " << crossProdVec[3] << " " << crossProdVec[4]
                  << " cor(Au,u)=" << arma::as_scalar(arma::cor(crossProdVec, bVec))
                  << std::endl;
        std::string dump_path = "/tmp/cpp_Ku_call" + std::to_string(crossprod_call_count) + ".csv";
        std::ofstream df(dump_path);
        if (df.is_open()) {
            df << "u,Au\n";
            for (arma::uword k = 0; k < crossProdVec.n_elem; ++k)
                df << bVec[k] << "," << crossProdVec[k] << "\n";
            df.close();
            std::cout << "  [C++ dumped " << crossProdVec.n_elem << " rows to " << dump_path << "]" << std::endl;
        }
        // Dump ptrsubSampleInGeno: FAM row index (1-based) for each GRM sample index
        if (crossprod_call_count == 0) {
            std::ofstream pf("/tmp/cpp_ptrsub.csv");
            if (pf.is_open()) {
                pf << "grm_idx,fam_row_1based\n";
                for (size_t i = 0; i < geno.ptrsubSampleInGeno.size(); i++)
                    pf << i << "," << geno.ptrsubSampleInGeno[i] << "\n";
                pf.close();
                std::cout << "  [C++ dumped ptrsubSampleInGeno (" << geno.ptrsubSampleInGeno.size()
                          << " entries) first 10 fam rows: ";
                for (int i = 0; i < 10 && i < (int)geno.ptrsubSampleInGeno.size(); i++)
                    std::cout << geno.ptrsubSampleInGeno[i] << " ";
                std::cout << std::endl;
            }
        }
        // Per-marker dump: recompute vec + val1 for markers 0..2 serially
        {
            std::string mpath = "/tmp/cpp_marker_call" + std::to_string(crossprod_call_count) + ".csv";
            std::ofstream mf(mpath);
            if (mf.is_open()) {
                mf << "marker,sample,vec,u,val1\n";
                arma::fvec vec;
                for (int m = 0; m < 3; m++) {
                    geno.Get_OneSNP_StdGeno(m, &vec);
                    float val1 = arma::dot(vec, bVec);
                    std::cout << "  [C++ marker " << m << "] freq=" << geno.alleleFreqVec[m]
                              << " invstd=" << geno.invstdvVec[m]
                              << " |vec|=" << arma::norm(vec)
                              << " val1=dot(vec,u)=" << val1
                              << " vec[0:5]=" << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << " " << vec[4]
                              << std::endl;
                    for (arma::uword s = 0; s < vec.n_elem; s++)
                        mf << m << "," << s << "," << vec[s] << "," << bVec[s] << "," << val1 << "\n";
                }
                mf.close();
            }
        }
        crossprod_call_count++;
    }
#else
    if (crossprod_call_count < 3) crossprod_call_count++;
#endif

  	return(crossProdVec);
}


// Phase-2: multi-RHS analogue of getCrossprodMatAndKin. ψ·B for N×k B.
// Sparse-kinship branch does one sparse×dense product; dense branch uses the
// batched TBB/AVX2 path.
arma::fmat getCrossprodMatAndKinMat(const arma::fmat& Bmat){
	if (isUseSparseSigmaforInitTau | isUseSparseSigmaforModelFitting) {
		arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
		arma::mat x = result * arma::conv_to<arma::mat>::from(Bmat);
		return arma::conv_to<arma::fmat>::from(x);
	}
	return parallelCrossProdMat(Bmat);
}


// Phase-2: multi-RHS analogue of getCrossprod. Σ·P for N×k P where
// Σ = tau0·diag(1/w) + tau1·ψ.
arma::fmat getCrossprodMat(const arma::fmat& Pmat, arma::fvec& wVec,
                           arma::fvec& tauVec){
	if (tauVec(1) == 0) {
		arma::fmat out = Pmat;
		out.each_col() %= (tauVec(0) / wVec);
		return out;
	}
	arma::fmat crossProd1 = getCrossprodMatAndKinMat(Pmat);
	arma::fmat out = Pmat;
	out.each_col() %= (tauVec(0) / wVec);
	out += tauVec(1) * crossProd1;
	return out;
}


// Phase-2 rollback switch: SAIGE_NO_BLOCKPCG=1 restores the serial per-probe
// (GetTrace/GetTrace_q) and per-marker (VR) PCG loops.
bool isBlockPCGdisabled() {
	static const bool v = [](){
		const char* e = std::getenv("SAIGE_NO_BLOCKPCG");
		return e && std::string(e) == "1"; }();
	return v;
}


// INTERNAL: LOCO version of cross product with kinship matrix
arma::fvec getCrossprodMatAndKin_LOCO(arma::fcolvec& bVec){

        arma::fvec crossProdVec = parallelCrossProd_LOCO(bVec) ;
        //arma::fvec crossProdVec_2 = parallelCrossProd_LOCO_2(bVec) ;

	//for(int k=0; k < 10; k++) {
        //	std::cout << "new crossProdVec " << k << " " << crossProdVec[k] << std::endl;
        //	std::cout << "old crossProdVec " << k << " " << crossProdVec_2[k] << std::endl;

	//}	


        return(crossProdVec);
}


// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
struct indicesRelatedSamples : public RcppParallel::Worker {

  int  Ntotal;
  std::vector< std::pair<int, int> > &output;
  std::vector<float> &kinValues;  // FIX: Also store kinship values
  std::mutex output_mutex;

  indicesRelatedSamples(int Ntotal, std::vector< std::pair<int, int> > &output, std::vector<float> &kinValues) :
    Ntotal(Ntotal), output(output), kinValues(kinValues) {}


  void operator()(std::size_t begin, size_t end) {
    int m_M_Submarker = getSubMarkerNum();
    for(std::size_t k=begin; k < end; k++) {
      int i = (int)(k / Ntotal);
      int j = (int)(k % Ntotal);
      if((j <= i)){
                        i = Ntotal - i - 2;
                        j = Ntotal - j - 1;
      }
      //std::cout << "i,j,k debug: " << i << " " << j << " " << k << std::endl;
      // DEBUG: Check matrix column bounds before arma::dot
      if(i >= (int)geno.stdGenoMultiMarkersMat.n_cols || j >= (int)geno.stdGenoMultiMarkersMat.n_cols) {
          std::cerr << "[ERROR findIndiceRelatedSample] i=" << i << " j=" << j
                    << " n_cols=" << geno.stdGenoMultiMarkersMat.n_cols
                    << " n_rows=" << geno.stdGenoMultiMarkersMat.n_rows
                    << " k=" << k << " Ntotal=" << Ntotal << std::endl;
      }
      if(i < 0 || j < 0) {
          std::cerr << "[ERROR findIndiceRelatedSample] NEGATIVE INDEX i=" << i << " j=" << j
                    << " k=" << k << " Ntotal=" << Ntotal << std::endl;
      }
      float kinValueTemp = arma::dot((geno.stdGenoMultiMarkersMat).col(i), (geno.stdGenoMultiMarkersMat).col(j));
      kinValueTemp = kinValueTemp/m_M_Submarker;
      if(kinValueTemp >=  geno.relatednessCutoff) {
        std::lock_guard<std::mutex> lock(output_mutex);
        output.push_back( std::pair<int, int>(i, j) );
        kinValues.push_back(kinValueTemp);  // FIX: Store kinship value
      }
    }
  }

};


// INTERNAL: Utility function for printing combination indices (debug)
void printComb(int N){
  int x = N*(N-1)/2 - 1;
  for(std::size_t k=0; k < x; k++) {
      int i = k / N;
      int j = k % N;
      if((j < i)){
                        i = N - i - 2;
                        j = N - j - 1;
      }
     std::cout << "i,j " << i << "," << j << std::endl;
  }

}


//arma::fmat findIndiceRelatedSample(){
//arma::fmat findIndiceRelatedSample(){

// INTERNAL: Utility function to identify related samples based on kinship threshold
void findIndiceRelatedSample(){

  int Ntotal = geno.getNnomissing();
//  tbb::concurrent_vector< std::pair<float, float> > output;

//  indicesRelatedSamples indicesRelatedSamples(Ntotal,output);
  geno.indiceVec.clear();  // Clear before populating
  geno.kinValueVecSparse.clear();  // Clear kinship values
  indicesRelatedSamples indicesRelatedSamples(Ntotal, geno.indiceVec, geno.kinValueVecSparse);

  long int Ntotal2 = (long int)Ntotal;

  long int totalCombination = Ntotal2*(Ntotal2-1)/2 - 1;
  std::cout << "Ntotal: " << Ntotal << std::endl;
  std::cout << std::numeric_limits<int>::max() << std::endl;
  std::cout << std::numeric_limits<long int>::max() << std::endl;
  std::cout << std::numeric_limits<long long int>::max() << std::endl;
  std::cout << "totalCombination: " << totalCombination << std::endl;
  long int x = 1000001;
  int b = (int)(x / Ntotal);
  int a = (int)(x % Ntotal);
  std::cout << "a " << a << std::endl;
  std::cout << "b " << b << std::endl;
  
  parallelFor(0, totalCombination, indicesRelatedSamples);

//  arma::fmat xout(output.size()+Ntotal,2);

//  for(int i=0; i<output.size(); i++) {
//    xout(i,0) = output[i].first;
//    xout(i,1) = output[i].second;
//  }
//  for(int i=output.size(); i < output.size()+Ntotal; i++) {
//    xout(i,0) = i - output.size();
//    xout(i,1) = xout(i,0);
//  }

/*
  for(int i=0; i < Ntotal; i++){
    (geno.indiceVec).push_back( std::pair<int, int>(i, i) );
  }
*/

//  return(xout);
}



struct sparseGRMUsingOneMarker : public Worker {
   // input matrix to read from
  // arma::imat & iMat;
   // output matrix to write to
   arma::fvec & GRMvec;

   //int M = geno.getM();
   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
//   sparseGRMUsingOneMarker(arma::imat & iMat, arma::fvec &GRMvec)
//      : iMat(iMat), GRMvec(GRMvec) {}


  sparseGRMUsingOneMarker(arma::fvec &GRMvec)
      : GRMvec(GRMvec) {}


   // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
//            int iint = iMat(i,0);
//            int jint = iMat(i,1);
	   int iint = (geno.indiceVec)[i].first;	
	   int jint = (geno.indiceVec)[i].second;	
/*
            float ival = geno.m_OneSNP_StdGeno(iint);
            float jval = geno.m_OneSNP_StdGeno(jint);
            // write to output matrix
            //rmat(i,j) = sqrt(.5 * (d1 + d2));
            GRMvec(i) = ival*jval/M;
*/
	//use Look-Up table for calucate GRMvec(i)
	    int ival = geno.m_OneSNP_Geno(iint);	
	    int jval = geno.m_OneSNP_Geno(jint);
	    GRMvec(i) = geno.sKinLookUpArr[ival][jval]; 

      }
   }
};

//void parallelcalsparseGRM(arma::imat & iMat, arma::fvec &GRMvec) {

// INTERNAL: Parallel computation of sparse genetic relationship matrix
void parallelcalsparseGRM(arma::fvec &GRMvec) {

//  int n1 = geno.indiceVec.size();
  // allocate the output matrix
  //GRMvec.set_size(n1);
//  std::cout << "OKKK3: "  << std::endl;
//  sparseGRMUsingOneMarker sparseGRMUsingOneMarker(iMat, GRMvec);
  sparseGRMUsingOneMarker sparseGRMUsingOneMarker(GRMvec);
//  std::cout << "OKKK4: "  << std::endl;

//  std::cout << "n1 " << n1 << std::endl;
//  std::cout << "iMat.n_cols " << iMat.n_cols << std::endl;
  // call parallelFor to do the work
//  parallelFor(0, iMat.n_rows, sparseGRMUsingOneMarker);
  parallelFor(0, (geno.indiceVec).size(), sparseGRMUsingOneMarker);

  // return the output matrix
  // return GRMvec;
}


struct sumTwoVec : public Worker
{   
   // source vectors
   arma::fvec &x;
   
   arma::fvec &sumVec;
  
   //int M = geno.getM(); 
   // constructors
   sumTwoVec(arma::fvec &x,arma::fvec &sumVec) 
      : x(x), sumVec(sumVec) {}
   
     // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
            sumVec(i) = x(i)+(geno.kinValueVecFinal)[i];
	    (geno.kinValueVecFinal)[i] = sumVec(i);	
      }
   }
   
};

// INTERNAL: Parallel utility for summing two vectors  
void  parallelsumTwoVec(arma::fvec &x) {

  int n1 = x.n_elem;
  // allocate the output matrix
  arma::fvec sumVec;
  sumVec.set_size(n1);

  sumTwoVec sumTwoVec(x, sumVec);

  // call parallelFor to do the work
  parallelFor(0, x.n_elem, sumTwoVec);

}




// R CONNECTION: Core initialization function called from SAIGE_fitNULLGLMM_fast() in R
// Sets up global geno object with PLINK files and sample information
void setgeno(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool isDiagofKinSetAsOne)
{
	int start_s=clock();
        geno.setGenoObj(bedfile, bimfile, famfile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne);
	//geno.printAlleleFreqVec();
	//geno.printGenoVec();
	int stop_s=clock();
	cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
}





// R CONNECTION: Returns raw genotype data for a single SNP to R functions
// Used in association testing and quality control procedures
arma::ivec Get_OneSNP_Geno(int SNPIdx)
{

	arma::ivec temp = * geno.Get_OneSNP_Geno(SNPIdx);
	return(temp);

}



arma::ivec Get_OneSNP_Geno_forVarRatio(int SNPIdx)
{
       
        arma::ivec temp = * geno.Get_OneSNP_Geno_forVarRatio(SNPIdx);
        return(temp);

}



  
// R CONNECTION: Returns standardized genotype data for a single SNP to R functions
// Standardized genotypes are used in statistical computations and association tests
arma::fvec Get_OneSNP_StdGeno(int SNPIdx)
{

	arma::fvec temp; 
	geno.Get_OneSNP_StdGeno(SNPIdx, & temp);
//	for(int j = 0; j < 100; j++){
//                std::cout << "temp(j): " << j << " " << temp(j) << std::endl;

 //       }


	return(temp);

}
  
    
  

//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
// INTERNAL: Compute diagonal elements of sigma matrix
arma::fvec getDiagOfSigma(arma::fvec& wVec, arma::fvec& tauVec){
#ifdef SAIGE_DEBUG_IO
  fprintf(stderr, "[DBG1b] getDiagOfSigma enter\n"); fflush(stderr);
#endif
	int Nnomissing = geno.getNnomissing();
	//int M = geno.getM();
	int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
	//cout << "MminMAF=" << MminMAF << endl;
	//cout << "M=" << M << endl;
	arma::fvec diagVec(Nnomissing);
	float diagElement;
	float floatBuffer;
  	//float minvElement;
  
  	if(!(geno.setKinDiagtoOne)){ 
	  //diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /M + tauVec(0)/wVec;
	  diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /MminMAF + tauVec(0)/wVec;

	}else{
	  diagVec = tauVec(1) + tauVec(0)/wVec;
	}

	//std::cout << "M " << M << std::endl;
	//std::cout << "tauVec(0) " << tauVec(0) << std::endl;
	//std::cout << "tauVec(1) " << tauVec(1) << std::endl;
        //for(unsigned int i=0; i< 10; i++){
	//	 std::cout << "diagVec(i) " << diagVec(i) << std::endl;
	//}

	//make diag of kin to be 1 to compare results of emmax and gmmat
	//diagVec = tauVec(1) + tauVec(0)/wVec;
	for(unsigned int i=0; i< Nnomissing; i++){
//	if(i < 100){
//		std::cout << i << "th element of diag of sigma and wVec " << diagVec(i) << " " << wVec(i) << std::endl;
//	}
  		if(diagVec(i) < 1e-4){
  			diagVec(i) = 1e-4 ;
  		}
  	}
  


    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
  	return(diagVec);
}

// INTERNAL: LOCO version - compute diagonal elements of sigma matrix
arma::fvec getDiagOfSigma_LOCO(arma::fvec& wVec, arma::fvec& tauVec){

        int Nnomissing = geno.getNnomissing();
        int Msub = geno.getMsub();
        //cout << "N=" << N << endl;
        arma::fvec diagVec(Nnomissing);
        float diagElement;
        float floatBuffer;
        //float minvElement;
        diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO());
	int Msub_MAFge_minMAFtoConstructGRM_in_b = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
	int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = geno.getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
	
	diagVec = diagVec/(Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b) + tauVec(0)/wVec;
        //diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO()) /(Msub_MAFge_minMAFtoConstructGRM) + tauVec(0)/wVec;
        for(unsigned int i=0; i< Nnomissing; i++){
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
        return(diagVec);

}


// R CONNECTION: Returns diagonal elements of covariance matrix Sigma for survival analysis to R functions
// Used in survival GWAS for variance component estimation and statistical inference
arma::fvec getDiagOfSigma_surv(arma::fvec& diagofWminusUinv, arma::fvec& tauVec){

        int Nnomissing = geno.getNnomissing();
        int M = geno.getM();
        int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();

        //cout << "N=" << N << endl;
        arma::fvec diagVec(Nnomissing);
        float diagElement;
        float floatBuffer;
        //float minvElement;

        if(!(geno.setKinDiagtoOne)){
          diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /MminMAF + tauVec(0) * diagofWminusUinv;

        }else{
          diagVec = tauVec(1) + tauVec(0) * diagofWminusUinv;
        }
        //cout << "wVec: " << endl;
        //wVec.print();
        //std::cout << "M " << M << std::endl;
        //std::cout << "tauVec(0) " << tauVec(0) << std::endl;
        //std::cout << "tauVec(1) " << tauVec(1) << std::endl;
       // for(unsigned int i=0; i< Nnomissing; i++){
        //      std::cout << "(*geno.Get_Diagof_StdGeno()) /M: " << (*geno.Get_Diagof_StdGeno()) /M << std::endl;
        //}

        //make diag of kin to be 1 to compare results of emmax and gmmat
        //diagVec = tauVec(1) + tauVec(0)/wVec;
        for(unsigned int i=0; i< Nnomissing; i++){
//      if(i < 100){
//              std::cout << i << "th element of diag of sigma and wVec " << diagVec(i) << " " << wVec(i) << std::endl;
//      }
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }



    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
        return(diagVec);
}


// R CONNECTION: LOCO version of diagonal Sigma elements for survival analysis to R functions
// Used in leave-one-chromosome-out survival analysis to avoid genomic inflation
arma::fvec getDiagOfSigma_surv_LOCO(arma::fvec& diagofWminusUinv, arma::fvec& tauVec){

        int Nnomissing = geno.getNnomissing();
        int Msub = geno.getMsub();

        //cout << "N=" << N << endl;
        arma::fvec diagVec(Nnomissing);
        float diagElement;
        float floatBuffer;
        //float minvElement;
        int Msub_MAFge_minMAFtoConstructGRM = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
        diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO()) /(Msub_MAFge_minMAFtoConstructGRM) + tauVec(0)*diagofWminusUinv;
        for(unsigned int i=0; i< Nnomissing; i++){
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
        return(diagVec);

}




arma::fcolvec getCrossprod(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
	
	//for(int j = 0; j < 100; j++){
        //        std::cout << "bVec(j): " << bVec(j) << std::endl;
        //        std::cout << "crossProd1(j): " << crossProd1(j) << std::endl;

        //}


        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

	//for(int j = 0; j < 10; j++){
        //        std::cout << "crossProdVec(j): " << j << " " << crossProdVec(j) << std::endl;

        //}



        return(crossProdVec);
}




arma::fcolvec getCrossprod_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);
        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

        return(crossProdVec);
}


// INTERNAL: Extract vector elements at specific time point for survival analysis
arma::fvec extractVecatTimek(unsigned int ktime , arma::fvec & rvecIndex, arma::fvec & winvn) {
        arma::fvec kthVec;
        unsigned int n_kthVec = winvn.n_elem;
        kthVec.zeros(n_kthVec);
        for(unsigned int i=0; i< n_kthVec; i++){
                if(rvecIndex(i) >= ktime){
                        kthVec(i) = winvn(i);
                }
        }
        return(kthVec);
}

//extractUvecforkthTime(i, n_RvecIndex, n_NVec, n_sqrtDVec, vec);
// INTERNAL: Extract U vector for k-th time point in survival analysis
void extractUvecforkthTime(unsigned int kthtime, arma::fvec & RvecIndex,  arma::fvec& NVec,  arma::fvec & sqrtDVec, arma::fvec & kthVec){
        //unsigned int ktime=RvecIndex(nthsample);
        unsigned int nsample = RvecIndex.n_elem;
        kthVec.zeros();
        float sqrtDKth = sqrtDVec(kthtime);
        for(unsigned int j = 0; j < nsample; j++){
                if((RvecIndex(j)-1) >= kthtime){
                        kthVec(j) = NVec(j);
                }
        }
        kthVec = kthVec * sqrtDKth;
}


//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_UandbVec_surv : public Worker
{
        // source vectors
        arma::fcolvec n_bVec;
        arma::fcolvec n_RvecIndex;
        arma::fcolvec n_NVec;
        arma::fcolvec n_sqrtDVec;
        //unsigned int k_uniqTime;
        // product that I have accumulated
        arma::fvec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_UandbVec_surv(arma::fcolvec & x, arma::fvec & y,  arma::fvec & z,  arma::fvec & q)
                : n_bVec(x),n_RvecIndex(y),n_NVec(z),n_sqrtDVec(q) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(m_N);
        }
        CorssProd_UandbVec_surv(const CorssProd_UandbVec_surv& CorssProd_UandbVec_surv, Split)
                : n_bVec(CorssProd_UandbVec_surv.n_bVec),n_RvecIndex(CorssProd_UandbVec_surv.n_RvecIndex),n_NVec(CorssProd_UandbVec_surv.n_NVec),n_sqrtDVec(CorssProd_UandbVec_surv.n_sqrtDVec)
        {
                m_N = CorssProd_UandbVec_surv.m_N;
                m_bout.zeros(m_N);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                vec.zeros(m_N);
                //int nthsample;
                //int ktime;
                //arma::fvec vec.zeros(m_N);
                for(unsigned int i = begin; i < end; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractUvecforkthTime(i, n_RvecIndex, n_NVec, n_sqrtDVec, vec);
                        //for(unsigned int j = 0; j < ktime; j++){
                        //        vec(j) = n_Dvec(j)*n_sqrtWinvNVec(i);
                        //}
                        float val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_UandbVec_surv & rhs) {
                m_bout += rhs.m_bout;
        }
};




// R CONNECTION: Parallel computation of U matrix cross-product for survival analysis to R functions
// Optimized matrix operations using parallel processing for survival mixed models
arma::fvec parallelCrossProd_UandbVec_surv(arma::fcolvec & bVec, arma::fvec & RvecIndex, arma::fvec& NVec,  arma::fvec & sqrtDVec) {

//  // declare the InnerProduct instance that takes a pointer to the vector data
        unsigned int ktime = sqrtDVec.n_elem;
        CorssProd_UandbVec_surv  CorssProd_UandbVec_surv(bVec, RvecIndex, NVec, sqrtDVec);
        //int m_N = geno.getNnomissing();

//  // call paralleReduce to start the work
        parallelReduce(0, ktime, CorssProd_UandbVec_surv);

        return CorssProd_UandbVec_surv.m_bout;
}




// R CONNECTION: Computes (W-U)*b matrix product for survival analysis to R functions
// Essential operation in survival mixed model coefficient estimation and inference
arma::fcolvec getProdWminusUb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex, arma::fvec& NVec, arma::fvec& sqrtDVec, arma::fvec& wVec){
        //unsigned int nsample = geno.getNnomissing();
        //unsigned int kuniqtime = Dvec.n_elem;

        arma::fcolvec Ub = parallelCrossProd_UandbVec_surv(bVec, RvecIndex, NVec, sqrtDVec);
        arma::fcolvec WminusUb = wVec % bVec - Ub;
        return WminusUb;
}



// R CONNECTION: Computes cross-product operations for survival analysis covariance matrix to R functions
// Core computational function used in survival GWAS mixed model fitting
arma::fcolvec getCrossprod_Surv(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat & WinvNRt, arma::fmat & ACinv){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fmat WinvNRtG;
        arma::fmat ACivWinvNRtG;
        //cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
        //cout << "OKKKKK4" << endl;
        WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);

        crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        return(crossProdVec);
}



// R CONNECTION: LOCO version of survival cross-product operations to R functions
// Used in leave-one-chromosome-out survival analysis for unbiased association testing
arma::fcolvec getCrossprod_Surv_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat & WinvNRt, arma::fmat & ACinv){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fmat WinvNRtG;
        arma::fmat ACivWinvNRtG;
        //cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
        //cout << "OKKKKK4" << endl;
        WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);

        crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        return(crossProdVec);
}

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_WinvNRttandVec : public Worker
{
        // source vectors
        arma::fcolvec & n_bVec;
        arma::fcolvec & n_RvecIndex;
        arma::fcolvec & n_WinvN;

        unsigned int k_uniTime;

        // product that I have accumulated
        arma::fvec m_bout;


        // constructors
        CorssProd_WinvNRttandVec(arma::fcolvec & x, arma::fvec & y, arma::fvec & z, unsigned int k)
                : n_bVec(x),n_RvecIndex(y),n_WinvN(z),k_uniTime(k) {
                m_bout.zeros(k_uniTime);
        }
        CorssProd_WinvNRttandVec(const CorssProd_WinvNRttandVec& CorssProd_WinvNRttandVec, Split)
                : n_bVec(CorssProd_WinvNRttandVec.n_bVec),n_RvecIndex(CorssProd_WinvNRttandVec.n_RvecIndex),n_WinvN(CorssProd_WinvNRttandVec.n_WinvN),k_uniTime(CorssProd_WinvNRttandVec.k_uniTime)
        {

                m_bout.zeros(k_uniTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                float val1;
                int ktime;
                for(unsigned int i = begin; i < end; i++){
                        ktime = i;
                        vec=extractVecatTimek(ktime, n_RvecIndex, n_WinvN);
//                      std::cout << "j: " << j << std::endl;
                        val1 = dot(vec,  n_bVec);
                        m_bout[i] += m_bout[i] + val1;
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_WinvNRttandVec & rhs) {
        m_bout += rhs.m_bout;
        }
};

// INTERNAL: Extract vector for n-th sample in survival analysis
void extractVecfornthSample(unsigned int nthsample, unsigned int k_uniqTime, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & nthVec) {
        unsigned int ktime=RvecIndex(nthsample);
        nthVec.zeros(k_uniqTime);
        for(unsigned int j = 0; j < ktime; j++){
                nthVec(j) = sqrtWinvNVec(nthsample);
        }
}


// INTERNAL: Extract vector for n-th sample in survival analysis (double precision)
void extractVecfornthSample_double(unsigned int nthsample, unsigned int k_uniqTime, arma::vec & RvecIndex, arma::vec & sqrtWinvNVec, arma::vec & nthVec) {
        unsigned int ktime=RvecIndex(nthsample);
        nthVec.zeros(k_uniqTime);
        for(unsigned int j = 0; j < ktime; j++){
                nthVec(j) = sqrtWinvNVec(nthsample);
        }
}


//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_RandbVec_surv : public Worker
{
        // source vectors
        arma::fcolvec & n_bVec;
        arma::fcolvec & n_RvecIndex;
        unsigned int k_uniqTime;
        // product that I have accumulated
        arma::fvec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_RandbVec_surv(arma::fcolvec & x, arma::fvec & y, unsigned int k)
                : n_bVec(x),n_RvecIndex(y),k_uniqTime(k) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(k_uniqTime);
        }
        CorssProd_RandbVec_surv(const CorssProd_RandbVec_surv& CorssProd_RandbVec_surv, Split)
                : n_bVec(CorssProd_RandbVec_surv.n_bVec),n_RvecIndex(CorssProd_RandbVec_surv.n_RvecIndex),k_uniqTime(CorssProd_RandbVec_surv.k_uniqTime)
        {
                m_N = CorssProd_RandbVec_surv.m_N;
                m_bout.zeros(k_uniqTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                vec.zeros(k_uniqTime);
                float val1;
                int nthsample;
                int ktime;
                for(unsigned int i = begin; i < end; i++){
                        nthsample = i;
                        ktime=n_RvecIndex(i);
                        for(unsigned int j = 0; j < ktime; j++){
                                m_bout(j) += n_bVec(i);
                        }
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_RandbVec_surv & rhs) {
        m_bout += rhs.m_bout;
        }
};


//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_AandbVec_surv : public Worker
{
        // source vectors
        arma::fcolvec n_bVec;
        arma::fcolvec n_RvecIndex;
        arma::fcolvec n_sqrtWinvNVec;
        unsigned int k_uniqTime;
        // product that I have accumulated
        arma::fvec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_AandbVec_surv(arma::fcolvec & x, arma::fvec & y,  arma::fvec & z,  unsigned int k)
                : n_bVec(x),n_RvecIndex(y),n_sqrtWinvNVec(z),k_uniqTime(k) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(k_uniqTime);
        }
        CorssProd_AandbVec_surv(const CorssProd_AandbVec_surv& CorssProd_AandbVec_surv, Split)
                : n_bVec(CorssProd_AandbVec_surv.n_bVec),n_RvecIndex(CorssProd_AandbVec_surv.n_RvecIndex),n_sqrtWinvNVec(CorssProd_AandbVec_surv.n_sqrtWinvNVec),k_uniqTime(CorssProd_AandbVec_surv.k_uniqTime)
        {
                m_N = CorssProd_AandbVec_surv.m_N;
                m_bout.zeros(k_uniqTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                vec.zeros(k_uniqTime);
                //int nthsample;
                //int ktime;
                //arma::fvec vec.zeros(m_N);
                for(unsigned int i = begin; i < end; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractVecfornthSample(i, k_uniqTime, n_RvecIndex, n_sqrtWinvNVec, vec);
                        //for(unsigned int j = 0; j < ktime; j++){
                        //        vec(j) = n_Dvec(j)*n_sqrtWinvNVec(i);
                        //}
                        float val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_AandbVec_surv & rhs) {
                m_bout += rhs.m_bout;
        }
};

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_AandbVec_surv_double : public Worker
{
        // source vectors
        arma::colvec n_bVec;
        arma::colvec n_RvecIndex;
        arma::colvec n_sqrtWinvNVec;
        unsigned int k_uniqTime;
        // product that I have accumulated
        arma::vec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_AandbVec_surv_double(arma::colvec & x, arma::vec & y,  arma::vec & z,  unsigned int k)
                : n_bVec(x),n_RvecIndex(y),n_sqrtWinvNVec(z),k_uniqTime(k) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(k_uniqTime);
        }
        CorssProd_AandbVec_surv_double(const CorssProd_AandbVec_surv_double& CorssProd_AandbVec_surv_double, Split)
                : n_bVec(CorssProd_AandbVec_surv_double.n_bVec),n_RvecIndex(CorssProd_AandbVec_surv_double.n_RvecIndex),n_sqrtWinvNVec(CorssProd_AandbVec_surv_double.n_sqrtWinvNVec),k_uniqTime(CorssProd_AandbVec_surv_double.k_uniqTime)
        {
                m_N = CorssProd_AandbVec_surv_double.m_N;
                m_bout.zeros(k_uniqTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::vec vec;
                vec.zeros(k_uniqTime);
                //int nthsample;
                //int ktime;
                //arma::fvec vec.zeros(m_N);
                for(unsigned int i = begin; i < end; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractVecfornthSample_double(i, k_uniqTime, n_RvecIndex, n_sqrtWinvNVec, vec);
                        //for(unsigned int j = 0; j < ktime; j++){
                        //        vec(j) = n_Dvec(j)*n_sqrtWinvNVec(i);
                        //}
                        double val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_AandbVec_surv_double & rhs) {
                m_bout += rhs.m_bout;
        }
};



// R CONNECTION: Parallel computation of A matrix cross-product for survival analysis to R functions
// High-performance matrix operations for survival time-to-event modeling
arma::fvec parallelCrossProd_AandbVec_surv(arma::fcolvec & bVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, unsigned int kuniqtime) {

//  // declare the InnerProduct instance that takes a pointer to the vector data
//      unsigned int ktime = sqrtWinvNVec.n_elem;
        CorssProd_AandbVec_surv  CorssProd_AandbVec_surv(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
        int m_N = geno.getNnomissing();

//  // call paralleReduce to start the work
        parallelReduce(0, m_N, CorssProd_AandbVec_surv);

        return CorssProd_AandbVec_surv.m_bout;
}



arma::vec parallelCrossProd_AandbVec_surv_double(arma::colvec & bVec, arma::vec & RvecIndex, arma::vec & sqrtWinvNVec, unsigned int kuniqtime) {

//  // declare the InnerProduct instance that takes a pointer to the vector data
//      unsigned int ktime = sqrtWinvNVec.n_elem;
        CorssProd_AandbVec_surv_double  CorssProd_AandbVec_surv_double(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
        int m_N = geno.getNnomissing();

//  // call paralleReduce to start the work
        parallelReduce(0, m_N, CorssProd_AandbVec_surv_double);

        return CorssProd_AandbVec_surv_double.m_bout;
}




arma::fvec parallelCrossProd_RandbVec_surv(arma::fcolvec & bVec, arma::fvec & RvecIndex, unsigned int kuniqtime) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        CorssProd_RandbVec_surv  CorssProd_RandbVec_surv(bVec, RvecIndex, kuniqtime);
        int m_N = geno.getNnomissing();
  // call paralleReduce to start the work
        parallelReduce(0, m_N, CorssProd_RandbVec_surv);

        return CorssProd_RandbVec_surv.m_bout;
}


// R CONNECTION: Computes R*b matrix product for survival analysis to R functions
// Matrix operation for survival data risk set calculations and model fitting
arma::fcolvec getProdRb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex, unsigned int kuniqtime){
        //unsigned int nsample = geno.getNnomissing();
        //unsigned int kuniqtime = Dvec.n_elem;

        arma::fcolvec Rb = parallelCrossProd_RandbVec_surv(bVec, RvecIndex, kuniqtime);
        return Rb;
}



// R CONNECTION: Computes A*b matrix product for survival analysis to R functions
// Core matrix operation in survival mixed model variance component calculations
arma::fcolvec getProdAb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex, arma::fvec& sqrtWinvNVec,arma::fvec& Dvec){
        //unsigned int nsample = geno.getNnomissing();
        unsigned int kuniqtime = Dvec.n_elem;

        arma::fcolvec Ab = parallelCrossProd_AandbVec_surv(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
        //cout << "Ab 1st part " << endl;
        //Ab.print();
        Ab = Ab +  (-1/Dvec) % bVec;
        return Ab;
}


arma::colvec getProdAb_Surv_double(arma::colvec& bVec, arma::vec & RvecIndex, arma::vec& sqrtWinvNVec, arma::vec& Dvec){
        //unsigned int nsample = geno.getNnomissing();
        unsigned int kuniqtime = Dvec.n_elem;

        arma::colvec Ab = parallelCrossProd_AandbVec_surv_double(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
        Ab = Ab +  (-1/Dvec) % bVec;
        return Ab;
}



arma::fvec getDiagofA( arma::fvec& RvecIndex, arma::fvec& sqrtWinvNVec,arma::fvec& Dvec){
        arma::fvec diagA;
        diagA = (-1/Dvec);
        unsigned int nsample = sqrtWinvNVec.n_elem;
        arma::fvec vec;
        unsigned int k_uniqTime = Dvec.n_elem;

        for(unsigned int i = 0; i < nsample; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractVecfornthSample(i, k_uniqTime, RvecIndex, sqrtWinvNVec, vec);
                        //cout << "vec.n_elem: " << vec.n_elem << endl;
                        diagA = diagA + vec % vec;
        }


        for(unsigned int i = 0; i < k_uniqTime; i++){
                if(diagA(i) == 0.0){
                        diagA(i) = 0.0001;
                }
        }

        return(diagA);
}



arma::vec getDiagofA_double( arma::vec& RvecIndex, arma::vec& sqrtWinvNVec,arma::vec& Dvec){
        arma::vec diagA;
        diagA = (-1/Dvec);
        unsigned int nsample = sqrtWinvNVec.n_elem;
        arma::vec vec;
        unsigned int k_uniqTime = Dvec.n_elem;

        for(unsigned int i = 0; i < nsample; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractVecfornthSample_double(i, k_uniqTime, RvecIndex, sqrtWinvNVec, vec);
                        //cout << "vec.n_elem: " << vec.n_elem << endl;
                        diagA = diagA + vec % vec;
        }


        for(unsigned int i = 0; i < k_uniqTime; i++){
                if(diagA(i) == 0.0){
                        diagA(i) = 0.0001;
                }
        }

        return(diagA);

}



arma::vec getPCG1ofACinvAndVector_test(arma::vec& bVec,  arma::vec& RvecIndex, arma::vec& sqrtWinvNVec,arma::vec& Dvec, int maxiterPCG, float tolPCG, arma::vec & wVec, arma::vec & tauVec, arma::mat & Rmat){
    unsigned int kuniqtime = Dvec.n_elem;
    //cout << "kuniqtime is " << kuniqtime << endl;
    arma::vec xVec(kuniqtime);
    xVec.zeros();
        //bVec = bVec/(1e+3);
        arma::vec rVec = bVec;
        arma::vec r1Vec;
        arma::vec zVec(kuniqtime);
        arma::vec minvVec(kuniqtime);

        minvVec = 1/getDiagofA_double(RvecIndex,sqrtWinvNVec,Dvec);/////To update
        zVec = minvVec % rVec;
        cout << "minvVec(10): " << minvVec(10) << endl;
        cout << "minvVec(20): " << minvVec(20) << endl;
        //zVec = rVec;
        double sumr2 = sum(rVec % rVec);
        arma::vec z1Vec(kuniqtime);
        arma::vec pVec = zVec;
        cout << "Rmat.n_cols " << Rmat.n_cols << endl;
        cout << "Rmat.n_rows " << Rmat.n_rows << endl;
        cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
        //arma::fmat sqrtWinvNmat = arma::diagmat(sqrtWinvNVec);
        // cout << "OK" << endl;
        //arma::fmat ApVectemp = (Rmat.t()) * sqrtWinvNmat;
        //arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
        arma::colvec ApVec = getProdAb_Surv_double(pVec,RvecIndex,sqrtWinvNVec,Dvec);
        //cout << "ApVec(10): " << ApVec(10) << endl;
        //cout << "ApVec0(10): " << ApVec0(10) << endl;
        int iter = 0;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::colvec ApVec = getProdAb_Surv_double(pVec,RvecIndex,sqrtWinvNVec,Dvec);
                cout << "iter: " << iter << endl;
                //arma::fcolvec ApVectemp = (Rmat.t()) * sqrtWinvNVec;
                //arma::fcolvec ApVec = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
                //cout << "Rmat.n_cols " << Rmat.n_cols << endl;
                //cout << "Rmat.n_rows " << Rmat.n_rows << endl;
                //cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
                //arma::fmat ApVectemp = (Rmat.t()) * sqrtWinvNmat;
                //arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
                cout << "ApVec(10): " << ApVec(10) << endl;
                cout << "pVec(10): " << pVec(10) << endl;
                //cout << "ApVec0(10): " << ApVec0(10) << endl;
                arma::vec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);
                cout << "rVec.t() * zVec " << rVec.t() * zVec << endl;
                cout << "pVec.t() * ApVec " << pVec.t() * ApVec << endl;
                float a = preA(0);
                cout << "a: " << a << endl;
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                arma::vec z1Vec = minvVec % r1Vec;
                //arma::fvec z1Vec = r1Vec;
                arma::vec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                double bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                cout << "bet: " << bet << endl;
                cout << "Prebet.n_elem: " << Prebet.n_elem << endl;
                cout << "z1Vec(10): " << z1Vec(10) << endl;
                cout << "r1Vec(10): " << r1Vec(10) << endl;


                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
                cout << "sumr2 is " << sumr2 << endl;
        }
        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;
        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
        //xVec = xVec *(1e+3);
        return(xVec);
}


arma::fvec getPCG1ofACinvAndVector(arma::fvec& bVec,  arma::fvec& RvecIndex, arma::fvec& sqrtWinvNVec,arma::fvec& Dvec, int maxiterPCG, float tolPCG, arma::fvec & wVec, arma::fvec & tauVec){
    maxiterPCG = 200;
    unsigned int kuniqtime = Dvec.n_elem;
    //cout << "kuniqtime is " << kuniqtime << endl;
    arma::fvec xVec(kuniqtime);
    xVec.zeros();

        //bVec = bVec /(1e+3);
        arma::fvec rVec = bVec;
        arma::fvec r1Vec;
        arma::fvec zVec(kuniqtime);
        arma::fvec minvVec(kuniqtime);

        minvVec = 1/getDiagofA(RvecIndex,sqrtWinvNVec,Dvec);/////To update
        zVec = minvVec % rVec;
        //cout << "minvVec(10): " << minvVec(10) << endl;
        //cout << "minvVec(20): " << minvVec(20) << endl;
        //zVec = rVec;
        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(kuniqtime);
        arma::fvec pVec = zVec;

        //arma::fcolvec ApVec = getProdAb_Surv(pVec,RvecIndex,sqrtWinvNVec,Dvec);


        //cout << "RmatIndex.n_cols " << RmatIndex.n_cols << endl;
        //cout << "RmatIndex.n_rows " << RmatIndex.n_rows << endl;
        //cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
        //arma::fcolvec ApVectemp = (Rmat.t()) * sqrtWinvNVec;
        //arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
        //cout << "ApVec(10): " << ApVec(10) << endl;
        //cout << "ApVec0(10): " << ApVec0(10) << endl;



        int iter = 0;
        arma::fcolvec ApVec;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                ApVec = getProdAb_Surv(pVec,RvecIndex,sqrtWinvNVec,Dvec);
                //cout << "iter: " << iter << endl;
                //for(size_t j=0; j< 10; j++){
                //      cout << "j: " << j << " ApVec(j) " << ApVec(j) << endl;
                //}



                //arma::fcolvec ApVectemp = (Rmat.t()) * sqrtWinvNVec;
                //arma::fcolvec ApVec = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
                //cout << "RmatIndex.n_cols " << RmatIndex.n_cols << endl;
                //cout << "RmatIndex.n_rows " << RmatIndex.n_rows << endl;
                //cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
                //arma::fcolvec ApVectemp = (RmatIndex.t()) * sqrtWinvNVec;
                //arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec + Dvec % pVec;
                //cout << "ApVec(10): " << ApVec(10) << endl;
                //cout << "pVec(10): " << pVec(10) << endl;
                //cout << "ApVec0(0): " << ApVec0(10) << endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);
                float a = preA(0);
                //cout << "a: " << a << endl;
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                z1Vec = minvVec % r1Vec;
                //arma::fvec z1Vec = r1Vec;
                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                //cout << "bet: " << bet << endl;
                //cout << "Prebet.n_elem: " << Prebet.n_elem << endl;
                //cout << "z1Vec(10): " << z1Vec(10) << endl;
                //cout << "r1Vec(10): " << r1Vec(10) << endl;
                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }
        //if (iter >= maxiterPCG){
        //        cout << "pcg did not converge. You may increase maxiter number." << endl;
        //}
        //cout << "sumr2 is " << sumr2 << endl;
        //cout << "iter from getPCG1ofACinvAndVector " << iter << endl;
        //xVec = xVec *(1e+3);
        return(xVec);
}



arma::fcolvec getProdRtb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex){
        unsigned int kuniqtime = bVec.n_elem;
        arma::fcolvec bsumVec;
        arma::fcolvec Rtbvec;
        unsigned int m_N = geno.getNnomissing();
        Rtbvec.zeros(m_N);
        bsumVec.zeros(kuniqtime);
        bsumVec(0) = bVec(0);
        for(unsigned int i = 1; i < kuniqtime; i++){
                bsumVec(i) = bsumVec(i-1) + bVec(i);
        }
        int ktime;
        for(unsigned int j = 0; j < m_N; j++){
                ktime = RvecIndex(j);
                Rtbvec(j) = bsumVec(ktime-1);
        }
        return(Rtbvec);
}


// R CONNECTION: Optimized survival cross-product computation with PCG integration to R functions
// Enhanced version using preconditioned conjugate gradient for improved computational efficiency
arma::fcolvec getCrossprod_Surv_new(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & NWinv, arma::fvec & Dvec, unsigned int kuniqtime, int maxiterPCG, float tolPCG){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
        //cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
        //cout << "OKKKKK4" << endl;
        //cout << "crossProdVec0(0) " << crossProdVec0(0) << endl;
        //WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "NWinv: " << endl;
        //for(size_t i=0; i< 10; i++){

          //      cout << NWinv(i) << " " << endl;
        //}

        //cout << "bVec: " << endl;
        //for(size_t i=0; i< 10; i++){

          //      cout << bVec(i) << " " << endl;
        //}


        arma::fcolvec NWinvbVec =  NWinv % bVec;
        //cout << "OKKKKK5" << endl;



        /*
        cout << NWinvbVec.n_elem << endl;
        cout << NWinvbVec(0) << endl;
        cout << Rmat.n_cols << endl;
        cout << Rmat.n_rows << endl;

        arma::fmat Rmatt = Rmat.t();
        cout << Rmatt.n_rows << endl;
        cout << Rmatt.n_cols << endl;

        cout << "NWinvbVec: " << endl;
        for(size_t i=0; i< 10; i++){

                cout << NWinvbVec(i) << " " << endl;
        }



//      arma::fcolvec RNWinvb0 = Rmatt * NWinvbVec;


//      cout << "RNWinvb(0) " << RNWinvb(0) << endl;
*/
        RNWinvb = getProdRb_Surv(NWinvbVec, RvecIndex, kuniqtime);
 //     RNWinvb0 = (Rmat.t()) * NWinvbVec;
//      arma::fcolvec RNWinvb1 =   Rmat.t() * (NWinv % bVec);

//      cout << "RNWinvb(0) " << RNWinvb(0) << endl;
//      //cout << "RNWinvb0(0) " << RNWinvb0(0) << endl;
//      cout << "RNWinvb1(0) " << RNWinvb1(0) << endl;

//      cout << "OKKKKK5" << endl;
        //arma::fcolvec RNWinvb;
        // cout << "OKKKKK6" << endl;
        //cout << "RNWinvb(0) is " << RNWinvb(0) << endl;
        //arma::fcolvec DRNWinvb = Dvec % RNWinvb;
        //cout << "DRNWinvb(0) is " << DRNWinvb(0) << endl;
        // cout << "OKKKKK7" << endl;
        arma::fcolvec AinvRNWinvb;
          //for(size_t i=0; i< 5; i++){
          //                cout << "i: " << i << " RNWinvb(i) " << RNWinvb(i) << endl;
        //                   }
/*
          for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " RvecIndex(i) " << RvecIndex(i) << endl;
                                  }

          for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " sqrtWinvNVec(i) " << sqrtWinvNVec(i) << endl;
                                  }


          for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " Dvec(i) " << Dvec(i) << endl;
                                  }
        */
        float pxnorm = arma::norm(RNWinvb);

        RNWinvb = RNWinvb/pxnorm;
        AinvRNWinvb = getPCG1ofACinvAndVector(RNWinvb, RvecIndex, sqrtWinvNVec, Dvec, maxiterPCG, tolPCG, wVec, tauVec);
        AinvRNWinvb = AinvRNWinvb * pxnorm;
        //for(size_t i=0; i< 5; i++){
        //      cout << "i: " << i << " AinvRNWinvb(i) " << AinvRNWinvb(i) << endl;
        //}
        //arma::fmat sqrtWinvNRtDt(wVec.n_elem, Dvec.n_elem);
        //arma::fmat Dmat = diagmat(-1/Dvec);
        //arma::fmat sqrtWinvNmat = diagmat(sqrtWinvNVec);
        //cout << "OKKKKK7c" << endl;
        //arma::fmat sqrtWinvNRt0 = sqrtWinvNmat * Rmat;
        //cout << "OKKKKK7d" << endl;
        //cout << "sqrtWinvNRt0.n_cols: " << sqrtWinvNRt0.n_cols << endl;
        //cout << "sqrtWinvNVec.n_elem: " << sqrtWinvNVec.n_elem << endl;
        //arma::fmat sqrtWinvNRt2 = (sqrtWinvNRt0.t()) * sqrtWinvNRt0;
        //cout << "OKKKKK7e" << endl;
        //arma::fmat A = sqrtWinvNRt2 + Dmat;
        //cout << "A(0,0) " << A(0,0) << endl;
        //cout << "A(1,1) " << A(1,1) << endl;
        //cout << "A(2,2) " << A(2,2) << endl;

        //arma::fcolvec Adiag = getDiagofA(RvecIndex,sqrtWinvNVec,Dvec);
        //cout << "Adiag(0) " << Adiag(0)<< endl;
        //cout << "Adiag(1) " << Adiag(1)<< endl;
        //cout << "Adiag(1) " << Adiag(1)<< endl;

        //arma::fvec AinvRNWinvb0 = solve(A, RNWinvb);
        //cout << "AinvRNWinvb0(0) is " << AinvRNWinvb0(0) << endl;
        //cout << "AinvRNWinvb(0) is " << AinvRNWinvb(0) << endl;
        //cout << "AinvDRNWinvb(0) is " << AinvDRNWinvb(0) << endl;
        //cout << "AinvDRNWinvb0(0) is " << AinvDRNWinvb0(0) << endl;
        //cout << "OKKKKK7b" << endl;

        //arma::fcolvec DAinvDRNWinvb = Dvec % AinvDRNWinvb0;
        // cout << "OKKKKK8" << endl;
        // cout << "DAinvDRNWinvb(0) is " << DAinvDRNWinvb(0) << endl;
        //arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb0, RvecIndex);
        arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb, RvecIndex);
        //cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
        crossProdVec1 = NWinv % RtAinvDRNWinvb;
        //cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
        //cout << "OKKKKK9" << endl;
        //RtAinvDRNWinvb = Rmat/(1e+10)  * AinvRNWinvb;
        //cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
        //crossProdVec1 = NWinv % (Rmat  * AinvRNWinvb);
        //cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
        //crossProdVec1 = getprodWinvNRttandVec(ACivWinvNRtG, RmatIndex, WinvN, kuniqtime);
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }else{
                arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);

                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

                return(crossProdVec);
        }
}



arma::fvec getPCG1ofWminusUAndVector(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & NVec, arma::fvec & sqrtDVec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){

                   //  Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = sqrtDVec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;
    //cout << "xVec: " << endl;
    //xVec.print();
   // arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
        //cout << "rVec: " << endl;
        //rVec.print();
   arma::fvec rVec = bVec;

        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        //double wall1 = get_wall_time();
        //double cpu1  = get_cpu_time();
        //minvVec = diagofWminusUinv;
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
        //cout << "rVec1: " << endl;
        //dofWminusU.print();
        minvVec = 1/dofWminusU;
        //cout << "minvVec: " << endl;
        //minvVec.print();
        zVec = minvVec % rVec;
        //cout << "rVec: " << endl;
        //rVec.print();
        //cout << "rVec2: " << endl;
                //zVec = rVec;
        //double wall2 = get_wall_time();
        //double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
  //      cout << "sumr2: " << sumr2 << endl;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
                //cout << "OKKKKKK" << endl;

                //arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
                //arma::fcolvec RWinN =  Rmat.t() * WinvN;
                //cout << "RWinN(0) is " << RWinN(0) << endl;


                //cout << "RWinNpVec(0) is " << RWinNpVec(0) << endl;
                arma::fcolvec ApVec = getProdWminusUb_Surv(pVec, RvecIndex, NVec, sqrtDVec, wVec);
                //arma::fcolvec ApVec = getCrossprod_Surv_new(pVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
                //cout << "ApVec is " << ApVec(0) << endl;
                //cout << "OKKKKKK2" << endl;
                /*
                arma::fcolvec ApVec0;
                arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
                WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
        */




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        //if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        //}else{
        //        z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        //}

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
//        cout << "iter from getPCG1ofSigmaAndVector_WminusU " << iter << endl;
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}



// R CONNECTION: Second-generation survival cross-product computation to R functions
// Advanced implementation with enhanced numerical stability for survival mixed models
arma::fcolvec getCrossprod_Surv_new2(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex,  arma::fvec & NVec, arma::fvec & sqrtDVec, arma::fcolvec & diagofWminusUinv, unsigned int kuniqtime, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
        //arma::fcolvec sqrtDVec = arma::sqrt(Dvec);
        arma::fcolvec x0Vec(bVec.n_elem);
        x0Vec.zeros();
        //cout << "OKKKKK2" << endl;
        crossProdVec0 = tauVec(0)*getPCG1ofWminusUAndVector(wVec,tauVec,bVec,RvecIndex,NVec,sqrtDVec,diagofWminusUinv,x0Vec,maxiterPCG, tolPCG, dofWminusU);

        //cout << "OKKKKK3" << endl;
        if(tauVec(1) == 0){

                return(crossProdVec0);
        }else{
                arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);

                crossProdVec = crossProdVec0 +  tauVec(1)*crossProd1;

                return(crossProdVec);
        }
}




// R CONNECTION: LOCO version of second-generation survival cross-product to R functions
// Advanced LOCO implementation for leave-one-chromosome-out survival analysis
arma::fcolvec getCrossprod_Surv_new2_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex,  arma::fvec & NVec, arma::fvec & sqrtDVec,  arma::fcolvec & diagofWminusUinv, unsigned int kuniqtime, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        //arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
        //cout << "OKKKKK3" << endl;
        //arma::fcolvec sqrtDVec = arma::sqrt(Dvec);
        arma::fcolvec x0Vec(bVec.n_elem);
        x0Vec.zeros();
        crossProdVec0 = tauVec(0)*getPCG1ofWminusUAndVector(wVec,tauVec,bVec,RvecIndex,NVec,sqrtDVec,diagofWminusUinv,x0Vec,maxiterPCG, tolPCG, dofWminusU);

        if(tauVec(1) == 0){

                return(crossProdVec0);
        }else{
                arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);

                crossProdVec = crossProdVec0 +  tauVec(1)*crossProd1;

                return(crossProdVec);
        }
}


// R CONNECTION: LOCO version of optimized survival cross-product computation to R functions
// Leave-one-chromosome-out implementation with PCG integration for survival analysis
arma::fcolvec getCrossprod_Surv_new_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & NWinv, arma::fvec & Dvec, unsigned int kuniqtime, int maxiterPCG, float tolPCG){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
        //cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
        //cout << "OKKKKK4" << endl;
        //cout << "crossProdVec0(0) " << crossProdVec0(0) << endl;
        //WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "NWinv: " << endl;
        //for(size_t i=0; i< 10; i++){

          //      cout << NWinv(i) << " " << endl;
        //}

        //cout << "bVec: " << endl;
        //for(size_t i=0; i< 10; i++){

          //      cout << bVec(i) << " " << endl;
        //}


        arma::fcolvec NWinvbVec =  NWinv % bVec;
        //cout << "OKKKKK5" << endl;



        /*
        cout << NWinvbVec.n_elem << endl;
        cout << NWinvbVec(0) << endl;
        cout << Rmat.n_cols << endl;
        cout << Rmat.n_rows << endl;

        arma::fmat Rmatt = Rmat.t();
        cout << Rmatt.n_rows << endl;
        cout << Rmatt.n_cols << endl;

        cout << "NWinvbVec: " << endl;
        for(size_t i=0; i< 10; i++){

                cout << NWinvbVec(i) << " " << endl;
        }



//      arma::fcolvec RNWinvb0 = Rmatt * NWinvbVec;


//      cout << "RNWinvb(0) " << RNWinvb(0) << endl;
*/
        RNWinvb = getProdRb_Surv(NWinvbVec, RvecIndex, kuniqtime);
 //     RNWinvb0 = (Rmat.t()) * NWinvbVec;
//      arma::fcolvec RNWinvb1 =   Rmat.t() * (NWinv % bVec);

//      cout << "RNWinvb(0) " << RNWinvb(0) << endl;
//      //cout << "RNWinvb0(0) " << RNWinvb0(0) << endl;
//      cout << "RNWinvb1(0) " << RNWinvb1(0) << endl;

//      cout << "OKKKKK5" << endl;
        //arma::fcolvec RNWinvb;
        // cout << "OKKKKK6" << endl;
        //cout << "RNWinvb(0) is " << RNWinvb(0) << endl;
        //arma::fcolvec DRNWinvb = Dvec % RNWinvb;
        //cout << "DRNWinvb(0) is " << DRNWinvb(0) << endl;
        // cout << "OKKKKK7" << endl;
        arma::fcolvec AinvRNWinvb;
        /*  for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " RNWinvb(i) " << RNWinvb(i) << endl;
                                  }

          for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " RvecIndex(i) " << RvecIndex(i) << endl;
                                  }

          for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " sqrtWinvNVec(i) " << sqrtWinvNVec(i) << endl;
                                  }


          for(size_t i=0; i< 10; i++){
                          cout << "i: " << i << " Dvec(i) " << Dvec(i) << endl;
                                  }
        */
        float pxnorm = arma::norm(RNWinvb);
        RNWinvb = RNWinvb/pxnorm;
        AinvRNWinvb = getPCG1ofACinvAndVector(RNWinvb, RvecIndex, sqrtWinvNVec, Dvec, maxiterPCG, tolPCG, wVec, tauVec);
        AinvRNWinvb = AinvRNWinvb * pxnorm;
        //for(size_t i=0; i< 10; i++){
        //      cout << "i: " << i << " AinvRNWinvb(i) " << AinvRNWinvb(i) << endl;
        //}
        //arma::fmat sqrtWinvNRtDt(wVec.n_elem, Dvec.n_elem);
        //arma::fmat Dmat = diagmat(-1/Dvec);
        //arma::fmat sqrtWinvNmat = diagmat(sqrtWinvNVec);
        //cout << "OKKKKK7c" << endl;
        //arma::fmat sqrtWinvNRt0 = sqrtWinvNmat * Rmat;
        //cout << "OKKKKK7d" << endl;
        //cout << "sqrtWinvNRt0.n_cols: " << sqrtWinvNRt0.n_cols << endl;
        //cout << "sqrtWinvNVec.n_elem: " << sqrtWinvNVec.n_elem << endl;
        //arma::fmat sqrtWinvNRt2 = (sqrtWinvNRt0.t()) * sqrtWinvNRt0;
        //cout << "OKKKKK7e" << endl;
        //arma::fmat A = sqrtWinvNRt2 + Dmat;
        //cout << "A(0,0) " << A(0,0) << endl;
        //cout << "A(1,1) " << A(1,1) << endl;
        //cout << "A(2,2) " << A(2,2) << endl;

        //arma::fcolvec Adiag = getDiagofA(RvecIndex,sqrtWinvNVec,Dvec);
        //cout << "Adiag(0) " << Adiag(0)<< endl;
        //cout << "Adiag(1) " << Adiag(1)<< endl;
        //cout << "Adiag(1) " << Adiag(1)<< endl;

        //arma::fvec AinvRNWinvb0 = solve(A, RNWinvb);
        //cout << "AinvRNWinvb0(0) is " << AinvRNWinvb0(0) << endl;
        //cout << "AinvRNWinvb(0) is " << AinvRNWinvb(0) << endl;
        //cout << "AinvDRNWinvb(0) is " << AinvDRNWinvb(0) << endl;
        //cout << "AinvDRNWinvb0(0) is " << AinvDRNWinvb0(0) << endl;
        //cout << "OKKKKK7b" << endl;

        //arma::fcolvec DAinvDRNWinvb = Dvec % AinvDRNWinvb0;
        // cout << "OKKKKK8" << endl;
        // cout << "DAinvDRNWinvb(0) is " << DAinvDRNWinvb(0) << endl;
        //arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb0, RvecIndex);
        arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb, RvecIndex);
        //cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
        crossProdVec1 = NWinv % RtAinvDRNWinvb;
        //cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
        //cout << "OKKKKK9" << endl;
        //RtAinvDRNWinvb = Rmat/(1e+10)  * AinvRNWinvb;
        //cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
        //crossProdVec1 = NWinv % (Rmat  * AinvRNWinvb);
        //cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
        //crossProdVec1 = getprodWinvNRttandVec(ACivWinvNRtG, RmatIndex, WinvN, kuniqtime);
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);

        crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        return(crossProdVec);
}

/*

arma::fcolvec getCrossprod_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);
        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

        return(crossProdVec);
}
*/



// REMOVED: get_wall_time() - use getTime() from src/UTIL.cpp instead



// INTERNAL: Generate sparse genetic relationship matrix
arma::sp_mat gen_sp_GRM() {
    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}



// R CONNECTION: Generates sparse covariance matrix Sigma to R functions
// Creates sparse representation of mixed model covariance structure for efficient computation
arma::sp_mat gen_sp_Sigma(arma::fvec& wVec,  arma::fvec& tauVec){
   arma::fvec dtVec = (1/wVec) * (tauVec(0));
//   dtVec.print();
   arma::vec valueVecNew = valueVec * tauVec(1);

   int nnonzero = valueVec.n_elem;
   for(size_t i=0; i< nnonzero; i++){
     if(locationMat(0,i) == locationMat(1,i)){
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
       valueVecNew(i) = valueVecNew(i) + dtVec(locationMat(0,i));
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
        if(valueVecNew(i) < 1e-4){
                        valueVecNew(i) = 1e-4 ;
                }


     }
   }

    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVecNew, dimNum, dimNum);
//    std::cout << "result.n_rows " << result.n_rows << std::endl;
//    std::cout << "result.n_cols " << result.n_cols << std::endl;
    //result.print();
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}




// REMOVED: get_cpu_time() - use getTime() from src/UTIL.cpp instead


/*

// R CONNECTION: Generates sparse genomic relationship matrix (GRM) to R functions
// Creates sparse kinship matrix for mixed model genetic relationship modeling
arma::sp_mat gen_sp_GRM() {
    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}

// INTERNAL: Generate sparse sigma matrix for computations
arma::sp_mat gen_sp_Sigma(arma::fvec& wVec,  arma::fvec& tauVec){
   arma::fvec dtVec = (1/wVec) * (tauVec(0));
//   dtVec.print();
   arma::vec valueVecNew = valueVec * tauVec(1);

   int nnonzero = valueVec.n_elem;
   for(size_t i=0; i< nnonzero; i++){
     if(locationMat(0,i) == locationMat(1,i)){
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
       valueVecNew(i) = valueVecNew(i) + dtVec(locationMat(0,i));
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
	if(valueVecNew(i) < 1e-4){
  			valueVecNew(i) = 1e-4 ;
  		}


     }
   }

    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVecNew, dimNum, dimNum);
//    std::cout << "result.n_rows " << result.n_rows << std::endl;
//    std::cout << "result.n_cols " << result.n_cols << std::endl;
    //result.print();
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}

*/


// R CONNECTION: Sparse linear system solver version 3 to R functions
// Solves sparse matrix systems using optimized algorithms for computational efficiency
arma::vec gen_spsolve_v3(arma::vec & yvec){
    // sparse x sparse -> sparse
    //arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    //arma::vec y = arma::linspace<arma::vec>(0, 5, dimNum);
    arma::sp_mat result = gen_sp_GRM();

    std::cout << "yvec.n_elem: " << yvec.n_elem << std::endl;
    std::cout << "result.n_rows: " << result.n_rows << std::endl;
    std::cout << "result.n_cols: " << result.n_cols << std::endl;
    arma::vec x = arma::spsolve(result, yvec);

    return x;
}


arma::fvec gen_spsolve_v4(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec & yvec){

    arma::vec yvec2 = arma::conv_to<arma::vec>::from(yvec);

    arma::sp_mat result = gen_sp_Sigma(wVec, tauVec);
#ifdef SAIGE_DEBUG_IO
    fprintf(stderr, "[DBG3] gen_sp_Sigma done nnz=%llu\n", (unsigned long long)result.n_nonzero); fflush(stderr);
#endif

    arma::vec x = arma::spsolve(result, yvec2);
#ifdef SAIGE_DEBUG_IO
    fprintf(stderr, "[DBG4] spsolve done\n"); fflush(stderr);
#endif

//double wall3in = get_wall_time();
// double cpu3in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;


    arma::fvec z = arma::conv_to<arma::fvec>::from(x);

//double wall4in = get_wall_time();
// double cpu4in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall4in - wall3in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu4in - cpu3in  << endl;


    return z;
}


//bool isUsePrecondM = false;
//bool isUseSparseSigmaforInitTau = false;



// INTERNAL: Set flag for using preconditioned matrix
void setisUsePrecondM(bool isUseSparseSigmaforPCG){
	isUsePrecondM = isUseSparseSigmaforPCG;
}

// INTERNAL: Set flag for using sparse sigma in initial tau estimation
void setisUseSparseSigmaforInitTau(bool isUseSparseSigmaforInitTau0){
	isUseSparseSigmaforInitTau = isUseSparseSigmaforInitTau0;
}



// INTERNAL: Set flag for using sparse sigma in null model fitting
void setisUseSparseSigmaforNullModelFitting(bool isUseSparseSigmaforModelFitting0){
        isUseSparseSigmaforModelFitting = isUseSparseSigmaforModelFitting0;
}

// INTERNAL: Set flag for using PCG with sparse sigma
void setisUsePCGwithSparseSigma(bool isUsePCGwithSparseSigma0){
         isUsePCGwithSparseSigma = isUsePCGwithSparseSigma0;
}


//Modified on 11-28-2018 to allow for a preconditioner for CG (the sparse Sigma)                                                                                                                                     //Sigma = tau[1] * diag(1/W) + tau[2] * kins
//This function needs the function getDiagOfSigma and function getCrossprod


// R CONNECTION: Core PCG solver called from getCoefficients() and used throughout R functions
// Implements preconditioned conjugate gradient to solve Sigma^(-1) * b efficiently
arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& wVec,
                                   const arma::fvec& tauVec,
                                   const arma::fvec& bVec,
                                   int maxiterPCG, float tolPCG)
{
#ifdef SAIGE_DEBUG_IO
    { const char _m[] = "[DBG1] PCG ENTER\n"; write(2, _m, sizeof(_m)-1); }
    fprintf(stderr, "[DBG1] PCG ENTER sparse=%d pcg=%d n=%zu\n",
            (int)isUseSparseSigmaforModelFitting, (int)isUsePCGwithSparseSigma,
            (size_t)bVec.n_elem); fflush(stderr);
#endif
    const arma::uword n = bVec.n_elem;
    if (n == 0) throw std::invalid_argument("PCG: bVec is empty");
    if (wVec.n_elem != n) throw std::invalid_argument("PCG: wVec.len != bVec.len");
    if (tauVec.n_elem < 2) throw std::invalid_argument("PCG: tauVec must have >=2 elems");

    // Direct sparse solve path (R default when usePCGwithSparseGRM=FALSE)
    // Matches R's getPCG1ofSigmaAndVector path 2:
    //   if (isUseSparseSigmaforModelFitting && !isUsePCGwithSparseSigma) → direct solve
    if (isUseSparseSigmaforModelFitting && !isUsePCGwithSparseSigma) {
        arma::fvec w   = wVec;
        arma::fvec tau = tauVec;
        arma::fvec b   = bVec;
        return gen_spsolve_v4(w, tau, b);
    }

    // make non-const copies to satisfy legacy APIs
    arma::fvec w  = wVec;
    arma::fvec tau = tauVec;

    arma::fvec xVec(n, arma::fill::zeros);
    arma::fvec rVec = bVec;
    arma::fvec zVec(n, arma::fill::zeros);
    arma::fvec minvVec(n, arma::fill::zeros);


    if (!isUsePrecondM) {
        // legacy takes non-const refs
        minvVec = 1.0f / getDiagOfSigma(w, tau);

        zVec    = minvVec % rVec;
    } else {
        zVec = gen_spsolve_v4(w, tau, rVec);
        if (zVec.n_elem != n)
            throw std::runtime_error("PCG: gen_spsolve_v4 returned wrong length");
    }


    arma::fvec pVec = zVec;
    float sumr2 = arma::dot(rVec, rVec);
    int   iter  = 0;


    while (sumr2 > tolPCG && iter < maxiterPCG) {
        ++iter;
        // getCrossprod expects non-const refs too
        arma::fcolvec ApVec = getCrossprod(pVec, w, tau);
        float a = arma::as_scalar((rVec.t() * zVec) / (pVec.t() * ApVec));

        xVec += a * pVec;

        arma::fvec r1Vec = rVec - a * ApVec;

        arma::fvec z1Vec;
        if (!isUsePrecondM) {
            z1Vec = minvVec % r1Vec;
        } else {
            z1Vec = gen_spsolve_v4(w, tau, r1Vec);
            if (z1Vec.n_elem != n)
                throw std::runtime_error("PCG: gen_spsolve_v4 (z1) wrong length");
        }

        float beta = arma::as_scalar((z1Vec.t() * r1Vec) / (zVec.t() * rVec));
        pVec = z1Vec + beta * pVec;
        zVec = std::move(z1Vec);
        rVec = std::move(r1Vec);
        sumr2 = arma::dot(rVec, rVec);
    }

    if (iter >= maxiterPCG)
        std::cout << "pcg did not converge (iter=" << iter << ")\n";
    else
        std::cout << "iter from getPCG1ofSigmaAndVector " << iter << "\n";

    return xVec;
}


// Phase-2 batched-RHS PCG: solves Sigma X = B for all columns of B at once.
// Structure mirrors the GPU branch's getPCGofSigmaAndMatrix
// (SAIGE-work/src/SAIGE_fitGLMM_fast.cpp): every column keeps its own
// alpha/beta/convergence state (batch of independent CGs, not block-CG), so
// per column the math is identical to getPCG1ofSigmaAndVector up to fp
// association order in ψ·B. Converged columns freeze (skip updates) but stay
// in the batch — the ψ·B cost at k≤64 is dominated by streaming the packed
// matrix, which shrinking the batch would not reduce.
arma::fmat getPCGofSigmaAndMatrix(const arma::fvec& wVec,
                                  const arma::fvec& tauVec,
                                  const arma::fmat& Bmat,
                                  int maxiterPCG, float tolPCG)
{
    const arma::uword N = Bmat.n_rows;
    const arma::uword k = Bmat.n_cols;

    // Paths without a batched implementation fall back per column: the
    // sparse direct solve, and the sparse preconditioner (isUsePrecondM —
    // silently switching it to Jacobi would change iteration counts and can
    // hit maxiterPCG).
    if (isUseSparseSigmaforModelFitting || isUsePrecondM) {
        arma::fmat X(N, k);
        for (arma::uword j = 0; j < k; ++j) {
            arma::fvec b = Bmat.col(j);
            X.col(j) = getPCG1ofSigmaAndVector(wVec, tauVec, b,
                                               maxiterPCG, tolPCG);
        }
        return X;
    }

    arma::fvec w   = wVec;
    arma::fvec tau = tauVec;

    arma::fmat X(N, k, arma::fill::zeros);
    arma::fmat R = Bmat;
    arma::fvec minvVec = 1.0f / getDiagOfSigma(w, tau);
    arma::fmat Z = R.each_col() % minvVec;
    arma::fmat P = Z;

    arma::fvec rz(k), sumr2(k);
    for (arma::uword j = 0; j < k; ++j) {
        rz(j)    = arma::dot(R.col(j), Z.col(j));
        sumr2(j) = arma::dot(R.col(j), R.col(j));
    }
    arma::uvec active(k);
    for (arma::uword j = 0; j < k; ++j) active(j) = (sumr2(j) > tolPCG);

    int iter = 0;
    while (arma::any(active) && iter < maxiterPCG) {
        iter++;
        arma::fmat AP = getCrossprodMat(P, w, tau);
        for (arma::uword j = 0; j < k; ++j) {
            if (!active(j)) continue;
            const float pAp = arma::dot(P.col(j), AP.col(j));
            const float a   = rz(j) / pAp;
            X.col(j) += a * P.col(j);
            R.col(j) -= a * AP.col(j);
            Z.col(j)  = minvVec % R.col(j);
            const float rz_new = arma::dot(R.col(j), Z.col(j));
            const float bta    = rz_new / rz(j);
            P.col(j)  = Z.col(j) + bta * P.col(j);
            rz(j)     = rz_new;
            sumr2(j)  = arma::dot(R.col(j), R.col(j));
            if (sumr2(j) <= tolPCG) active(j) = 0;
        }
    }
    if (arma::any(active)) {
        std::cout << "batched PCG: " << arma::sum(active) << "/" << k
                  << " columns did not converge in " << maxiterPCG
                  << " iterations\n";
    }
    std::cout << "iter from getPCGofSigmaAndMatrix " << iter << " for " << k
              << " RHS\n";
    return X;
}


// R CONNECTION: PCG solver for Sigma^(-1)*b in survival analysis to R functions
// Preconditioned conjugate gradient algorithm for survival mixed model linear systems
arma::fvec getPCG1ofSigmaAndVector_Surv(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
        int Nnomissing = geno.getNnomissing();
        arma::fvec xVec(Nnomissing);
        //xVec.zeros();
        xVec = x0Vec;
        //arma::fvec rVec(Nnomissing);
//if(isUseSparseSigmaforInitTau){
//        cout << "use sparse kinship to estimate initial tau " <<  endl;
//        xVec = gen_spsolve_v4(wVec, tauVec, bVec);
        if(isUseSparseSigmaforModelFitting){
             cout << "use sparse kinship to estimate the variance ratio " << endl;
        }
//      rVec = bVec - getCrossprod_Surv_sparseGRM(xVec, wVec, tauVec, WinvNRt, ACinv);
//}else{
        //arma::fvec rVec = bVec;
        arma::fvec rVec = bVec - getCrossprod_Surv(xVec, wVec, tauVec, WinvNRt, ACinv);
//}
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);

       double wall1 = get_wall_time();
       double cpu1  = get_cpu_time();

        if (!isUsePrecondM){
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                zVec = minvVec % rVec;
                //zVec = rVec;
        }else{


                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }


        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
        //cout << "OKKKKKK2" << endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "sumr2: " << sumr2 << std::endl;
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
//} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}



//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
//This function needs the function getDiagOfSigma and function getCrossprod

// R CONNECTION: Legacy PCG solver for Sigma^(-1)*b operations to R functions
// Original implementation maintained for compatibility and benchmarking
arma::fvec getPCG1ofSigmaAndVector_old(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG){
	           //  Start Timers
//    double wall0 = get_wall_time();
//    double cpu0  = get_cpu_time();

//	 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	//cout << "HELLO: "  << endl;
	//cout << "HELL2: "  << endl;
  	arma::fvec rVec = bVec;
	//cout << "HELLOa: "  << endl;
  	arma::fvec r1Vec;
	//cout << "HELLOb: "  << endl;
  	int Nnomissing = geno.getNnomissing();
	//cout << "HELL1: "  << endl;

  	arma::fvec crossProdVec(Nnomissing);
	//cout << "HELL2: "  << endl;

  	arma::fvec minvVec = 1/getDiagOfSigma(wVec, tauVec);
//	cout << "HELL3: "  << endl;
//	for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
  	float sumr2 = sum(rVec % rVec);
/*
	if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
  	arma::fvec zVec = minvVec % rVec;
  	arma::fvec z1Vec;
 	arma::fvec pVec = zVec;
	/*
        if(bVec[0] == 1 && bVec[2] == 1){
	for(int i = 0; i < 10; i++){ 
		cout << "pVec[i]: " << i << " " << pVec[i] << endl;
 	}
	}
*/
  	arma::fvec xVec(Nnomissing);
  	xVec.zeros();
  
  	int iter = 0;
  	while (sumr2 > tolPCG && iter < maxiterPCG) {
    		iter = iter + 1;
    		arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
    		arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

    		float a = preA(0);
	
/*	     if(bVec[0] == 1 && bVec[2] == 1){
			cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
        		for(int i = 0; i < 10; i++){
			
                		cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                		cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                		cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                		cout << "rVec[i]: " << i << " " << rVec[i] << endl;
        		}
        	    }   
*/	
 
    		xVec = xVec + a * pVec;
/*
		if(bVec[0] == 1 && bVec[2] == 1){
        		for(int i = 0; i < 10; i++){
                		cout << "xVec[i]: " << i << " " << xVec[i] << endl;
        		}
        	}   

*/

 
    		r1Vec = rVec - a * ApVec;
/*
		if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
    		z1Vec = minvVec % r1Vec;
    		arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
    		float bet = Prebet(0);
    		pVec = z1Vec+ bet*pVec;
    		zVec = z1Vec;
    		rVec = r1Vec;
    
    		sumr2 = sum(rVec % rVec);
/*
		if(bVec[0] == 1 && bVec[2] == 1){
			std::cout << "sumr2: " << sumr2 << std::endl;
			std::cout << "tolPCG: " << tolPCG << std::endl;
		}
*/
  	}
  
  	if (iter >= maxiterPCG){
    		cout << "pcg did not converge. You may increase maxiter number." << endl;
     
  	}
  	cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;

//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;
	
//	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
  	return(xVec);
}




// R CONNECTION: LOCO PCG solver for survival analysis to R functions
// Leave-one-chromosome-out preconditioned conjugate gradient for survival models
arma::fvec getPCG1ofSigmaAndVector_Surv_LOCO(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec,int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
        int Nnomissing = geno.getNnomissing();
        arma::fvec xVec(Nnomissing);
        //xVec.zeros();
        xVec = x0Vec;

if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else if(isUseSparseSigmaforModelFitting){
        cout << "use sparse kinship to fit the model " << endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else{
        arma::fvec rVec = bVec -  getCrossprod_Surv_LOCO(xVec, wVec, tauVec, WinvNRt, ACinv);
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);

       double wall1 = get_wall_time();
       double cpu1  = get_cpu_time();

        if (!isUsePrecondM){
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                minvVec = 1/getDiagOfSigma_surv_LOCO(diagofWminusUinv, tauVec);
                //zVec = minvVec % rVec;
                zVec = rVec;
        }else{


                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }


        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                arma::fcolvec ApVec = getCrossprod_Surv_LOCO(pVec, wVec, tauVec, WinvNRt, ACinv);
        //cout << "OKKKKKK2" << endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                //z1Vec = minvVec % r1Vec;
                z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "sumr2: " << sumr2 << std::endl;
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}



// R CONNECTION: Enhanced PCG solver for survival analysis to R functions
// Optimized preconditioned conjugate gradient with improved convergence
arma::fvec getPCG1ofSigmaAndVector_Surv_new(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = Dvec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;

    //if(isUseSparseSigmaforInitTau){
      //  cout << "use sparse kinship to estimate initial tau " <<  endl;
      //  xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
//    if(isUseSparseSigmaforModelFitting){
//      cout << "use sparse kinship to fit the model " << endl;
//      xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
//    }else{
        if(isUseSparseSigmaforModelFitting){
                cout << "use sparse kinship to estimate the variance ratio " << endl;
        }

        arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);

        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
                minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                zVec = minvVec % rVec;
                //zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
        double wall2 = get_wall_time();
        double cpu2  = get_cpu_time();
        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;

        int iter = 0;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;

                arma::fcolvec ApVec = getCrossprod_Surv_new(pVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);

                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

                xVec = xVec + a * pVec;


                r1Vec = rVec - a * ApVec;

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
        return(xVec);
}


// R CONNECTION: Second-generation PCG solver for survival analysis to R functions
// Advanced PCG implementation with enhanced numerical stability and performance
arma::fvec getPCG1ofSigmaAndVector_Surv_new2(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & NVec, arma::fvec & sqrtDvec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = sqrtDvec.n_elem;
    arma::fvec xVec(Nnomissing);
    xVec.zeros();
    //xVec = x0Vec;
    cout << "xVec: " << endl;
    //xVec.print();
    if(isUseSparseSigmaforModelFitting){
                 cout << "use sparse kinship to estimate the variance ratio " << endl;
     }

//    if(isUseSparseSigmaforInitTau){
//        cout << "use sparse kinship to estimate initial tau " <<  endl;
//        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
//    }else if(isUseSparseSigmaforModelFitting){
//      cout << "use sparse kinship to fit the model " << endl;
//      xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    //}else{
        //arma::fvec rVec = bVec;
        arma::fvec rVec = bVec - getCrossprod_Surv_new2(xVec,  wVec, tauVec, RvecIndex,NVec, sqrtDvec,diagofWminusUinv, kuniqtime, maxiterPCG, tolPCG, dofWminusU);
        //arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
        cout << "rVec: " << endl;
        //rVec.print();


        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
                minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                zVec = minvVec % rVec;
                //zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
        double wall2 = get_wall_time();
        double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
                //cout << "OKKKKKK" << endl;

                //arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
                //arma::fcolvec RWinN =  Rmat.t() * WinvN;
                //cout << "RWinN(0) is " << RWinN(0) << endl;


                //cout << "pVec(0) is " << pVec(0) << endl;
                arma::fcolvec ApVec = getCrossprod_Surv_new2(pVec, wVec, tauVec, RvecIndex,NVec, sqrtDvec, diagofWminusUinv, kuniqtime, maxiterPCG, tolPCG, dofWminusU);
                //cout << "ApVec is " << ApVec(0) << endl;


                //cout << "OKKKKKK2" << endl;
                /*
                arma::fcolvec ApVec0;
                arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
                WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
        */
                arma::fvec pAp = pVec.t() * ApVec;
                if(pAp(0) == 0){
                        return(xVec);
                }
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/

                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();
        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
       std::cout << "sumr2: " << sumr2 << std::endl;
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
//} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}


// R CONNECTION: LOCO second-generation PCG solver for survival analysis to R functions
// Advanced LOCO PCG with enhanced numerical stability for survival models
arma::fvec getPCG1ofSigmaAndVector_Surv_LOCO_new2(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & NVec, arma::fvec & sqrtDvec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = sqrtDvec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;
    //cout << "xVec: " << endl;
    //xVec.print();


    if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else if(isUseSparseSigmaforModelFitting){
        cout << "use sparse kinship to fit the model " << endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else{
        //arma::fvec rVec = bVec - getCrossprod_Surv_new2(xVec,  wVec, tauVec, RvecIndex,NVec, Dvec, kuniqtime, maxiterPCG, tolPCG);
        arma::fvec rVec = bVec - getCrossprod_Surv_new2_LOCO(xVec,  wVec, tauVec, RvecIndex,NVec, sqrtDvec, diagofWminusUinv, kuniqtime, maxiterPCG, tolPCG, dofWminusU);
        //arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
        //cout << "rVec: " << endl;
        //rVec.print();


        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
                minvVec = 1/getDiagOfSigma_surv_LOCO(diagofWminusUinv, tauVec);
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                zVec = minvVec % rVec;
                //zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
        double wall2 = get_wall_time();
        double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();
        int iter = 0;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
                //cout << "OKKKKKK" << endl;

                //arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
                //arma::fcolvec RWinN =  Rmat.t() * WinvN;
                //cout << "RWinN(0) is " << RWinN(0) << endl;


                //cout << "RWinNpVec(0) is " << RWinNpVec(0) << endl;
                arma::fcolvec ApVec = getCrossprod_Surv_new2_LOCO(pVec, wVec, tauVec, RvecIndex,NVec, sqrtDvec, diagofWminusUinv,  kuniqtime, maxiterPCG, tolPCG, dofWminusU);
                //cout << "ApVec is " << ApVec(0) << endl;
                //cout << "OKKKKKK2" << endl;
                /*
                arma::fcolvec ApVec0;
                arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
                WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
        */




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();
        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}



// R CONNECTION: Enhanced LOCO PCG solver for survival analysis to R functions
// Optimized leave-one-chromosome-out PCG for survival mixed models
arma::fvec getPCG1ofSigmaAndVector_Surv_new_LOCO(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec,  arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = Dvec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;

    if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else if(isUseSparseSigmaforModelFitting){
        cout << "use sparse kinship to fit the model " << endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else{
        arma::fvec rVec = bVec - getCrossprod_Surv_new_LOCO(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
                //minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                minvVec = 1/getDiagOfSigma_surv_LOCO(diagofWminusUinv, tauVec);
                zVec = minvVec % rVec;
                //zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
        double wall2 = get_wall_time();
        double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
        //cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
                //cout << "OKKKKKK" << endl;

                //arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
                //arma::fcolvec RWinN =  Rmat.t() * WinvN;
                //cout << "RWinN(0) is " << RWinN(0) << endl;

                arma::fcolvec ApVec = getCrossprod_Surv_new_LOCO(pVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
                //cout << "ApVec is " << ApVec(0) << endl;
                //cout << "OKKKKKK2" << endl;
                /*
                arma::fcolvec ApVec0;
                arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
                WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
        */




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}

//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
//This function needs the function getDiagOfSigma and function getCrossprod

// R CONNECTION: Standard LOCO PCG solver for mixed models to R functions
// Leave-one-chromosome-out preconditioned conjugate gradient for GLMM
arma::fvec getPCG1ofSigmaAndVector_LOCO(const arma::fvec& wVec_in,
                                        const arma::fvec& tauVec_in,
                                        const arma::fvec& bVec_in,
                                        int maxiterPCG, float tolPCG)
{
    // ---- shape checks ----
    const arma::uword n = bVec_in.n_elem;
    if (n == 0) throw std::invalid_argument("PCG_LOCO: bVec is empty");
    if (wVec_in.n_elem != n)
        throw std::invalid_argument("PCG_LOCO: wVec length (" + std::to_string(wVec_in.n_elem) +
                                    ") != bVec length (" + std::to_string(n) + ")");
    if (tauVec_in.n_elem < 2)
        throw std::invalid_argument("PCG_LOCO: tauVec must have at least 2 elements");

    // ---- make non-const local copies for legacy helpers ----
    arma::fvec wVec  = wVec_in;     // legacy APIs take arma::fvec&
    arma::fvec tauVec = tauVec_in;  // "
    arma::fvec bVec   = bVec_in;    // "

    // ---- init ----
    arma::fvec xVec(n, arma::fill::zeros);
    arma::fvec rVec = bVec;                 // residual
    arma::fvec zVec(n, arma::fill::zeros);
    arma::fvec minvVec(n, arma::fill::zeros);

    // LOCO behavior (which markers/blocks are excluded) should already be
    // configured globally via setStartEndIndex*()/set_Diagof_StdGeno_LOCO()
    // and is honored by the legacy helpers below.

    // Preconditioner.
    // R (src/SAIGE_fitGLMM_fast.cpp::getPCG1ofSigmaAndVector_LOCO) has NO sparse
    // branch here: it always uses the diagonal preconditioner built from the
    // LOCO diagonal. Mirror that exactly. (LOCO is disabled whenever the sparse
    // GRM is used to fit the null model, so the sparse branch is unreachable.)
    //
    // BUG FIX: this used to call the non-LOCO getDiagOfSigma()/getCrossprod(),
    // which made the "LOCO" solve numerically identical to the full-genome
    // solve — i.e. LOCO silently did nothing.
    minvVec = 1.0f / getDiagOfSigma_LOCO(wVec, tauVec);
    if (minvVec.n_elem != n)
        throw std::runtime_error("PCG_LOCO: getDiagOfSigma_LOCO returned wrong length");
    zVec = minvVec % rVec;

    arma::fvec pVec = zVec;
    float sumr2 = arma::dot(rVec, rVec);
    int   iter  = 0;

    while (sumr2 > tolPCG && iter < maxiterPCG) {
        ++iter;

        // Ap = Sigma_LOCO * p (the excluded chromosome block is subtracted
        // inside parallelCrossProd_LOCO)
        arma::fcolvec ApVec = getCrossprod_LOCO(pVec, wVec, tauVec);
        if (ApVec.n_elem != n)
            throw std::runtime_error("PCG_LOCO: getCrossprod_LOCO returned wrong length");

        float a = arma::as_scalar((rVec.t() * zVec) / (pVec.t() * ApVec));
        xVec += a * pVec;

        arma::fvec r1Vec = rVec - a * ApVec;

        arma::fvec z1Vec = minvVec % r1Vec;

        float beta = arma::as_scalar((z1Vec.t() * r1Vec) / (zVec.t() * rVec));
        pVec = z1Vec + beta * pVec;
        zVec = std::move(z1Vec);
        rVec = std::move(r1Vec);
        sumr2 = arma::dot(rVec, rVec);
    }

    if (iter >= maxiterPCG)
        std::cout << "pcg_loco did not converge (iter=" << iter << ")\n";
    else
        std::cout << "iter from getPCG1ofSigmaAndVector_LOCO " << iter << "\n";

    return xVec;   // length n
}


//http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/

// REMOVED: set_seed() - this function is commented out in src/UTIL.hpp and not implemented in src/UTIL.cpp

// REMOVED: nb() - use nb() from src/UTIL.cpp instead (takes unsigned int parameter)

// REMOVED: setChromosomeIndicesforLOCO() - the definition here was commented out
// while the header still declared it, so any caller would have failed to link.
// Restore definition + declaration together when LOCO is implemented.

// INTERNAL: Set start and end indices for chromosome analysis
void setStartEndIndex(int startIndex, int endIndex, int chromIndex){
  geno.startIndex = startIndex;
  geno.endIndex = endIndex;
  geno.Msub = 0;
  geno.chromIndex = chromIndex;

  for(size_t i=0; i< geno.M; i++){
	if(i < startIndex || i > endIndex){
  		if(geno.alleleFreqVec[i] >= minMAFtoConstructGRM && geno.alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
      
			geno.Msub = geno.Msub + 1;
  		}
	}
  }
  //geno.Msub = geno.M - (endIndex - startIndex + 1);
}



// INTERNAL: Set start and end index vectors for batch processing
void setStartEndIndexVec( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec){	
  geno.startIndexVec = startIndex_vec;
  geno.endIndexVec = endIndex_vec;
  //geno.Msub = geno.M - (endIndex - startIndex + 1);
}

// 
//void setStartEndIndexVec_forvr( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec){
//  geno.startIndexVec_forvr = startIndex_vec;
//  geno.endIndexVec_forvr = endIndex_vec;
  //geno.Msub = geno.M - (endIndex - startIndex + 1);
//}



//This function calculates the coefficients of variation for mean of a vector
// INTERNAL: Calculate coefficient of variation
float calCV(arma::fvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = arma::mean(xVec);
  float vecSd = arma::stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}

// Storage for pre-loaded random vectors from R (seed 200)
static std::vector<arma::fvec> preloaded_vectors;
static int preloaded_vector_idx = 0;
static bool use_preloaded_vectors = false;

// Load random vectors from CSV file
static void load_vectors_from_csv(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    std::cerr << "Warning: Could not open " << filepath << ", using random generation" << std::endl;
    use_preloaded_vectors = false;
    return;
  }

  preloaded_vectors.clear();
  std::string line;
  std::getline(file, line); // Skip header

  while (std::getline(file, line)) {
    std::vector<float> values;
    std::stringstream ss(line);
    std::string cell;
    std::getline(ss, cell, ','); // Skip vector_id
    while (std::getline(ss, cell, ',')) {
      values.push_back(std::stof(cell));
    }
    arma::fvec vec(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
      vec(i) = values[i];
    }
    preloaded_vectors.push_back(vec);
  }
  file.close();
  preloaded_vector_idx = 0;
  use_preloaded_vectors = true;
  std::cout << "Loaded " << preloaded_vectors.size() << " vectors from " << filepath << std::endl;
}

// Generate ±1 (Rademacher) vector of length n using C++ std library
// (arma::arma_rng::set_seed_random() can block waiting for entropy)
static inline arma::fvec rademacher_vec(int n) {
  // Check if we should use preloaded vectors
  if (use_preloaded_vectors && preloaded_vector_idx < (int)preloaded_vectors.size()) {
    arma::fvec vec = preloaded_vectors[preloaded_vector_idx];
    std::cout << "Using preloaded vector " << preloaded_vector_idx << std::endl;
    preloaded_vector_idx++;
    return vec;
  }

  // ERROR: If preloaded vectors are enabled but exhausted, stop!
  if (use_preloaded_vectors) {
    std::cerr << "ERROR: Ran out of preloaded vectors! Requested index=" << preloaded_vector_idx
              << ", but only " << preloaded_vectors.size() << " vectors loaded." << std::endl;
    std::cerr << "Set SAIGE_BYPASS_DIR to a directory with more precomputed random vectors." << std::endl;
    throw std::runtime_error("Preloaded vectors exhausted - cannot continue with matching R vectors");
  }

  // Fall back to random generation using R's RNG (matches R's Rcpp::rbinom(n,1,0.5))
  // Uses R C API directly to avoid Rcpp version mismatch in embedded mode
  // R's embedded runtime must be initialized (Rf_initEmbeddedR in main)
  // NOTE: GetRNGstate/PutRNGstate must NOT be called per-vector — R generates
  // all vectors in one continuous stream. The caller (GetTrace_q) calls set_seed()
  // which handles the RNG state. We just draw from the current state here.
  arma::fvec u(n);
  for (int i = 0; i < n; ++i) {
    double binom_val = Rf_rbinom(1.0, 0.5);
    u(i) = static_cast<float>(binom_val * 2.0 - 1.0);  // 0→-1, 1→+1 (matches R's uVec*2-1)
  }
  return u;
}

namespace saige {
float GetTrace(const arma::fmat& Sigma_iX,
               const arma::fmat& Xmat,
               const arma::fvec& wVec,
               const arma::fvec& tauVec,
               const arma::fmat& cov1,
               int nrun,
               int maxiterPCG,
               float tolPCG,
               float traceCVcutoff)
{
  std::cout << "=== Entering saige::GetTrace ===" << std::endl << std::flush;

  // Load precomputed vectors from R if file exists
  static bool vectors_loaded = false;
  static bool load_attempted = false;
  if (!load_attempted) {
    load_attempted = true;
    std::string bypass_path = saige_env_path("SAIGE_BYPASS_DIR", "random_vectors_seed10.csv");
    std::cout << "=== RANDOM VECTOR BYPASS CHECK ===" << std::endl;
    std::cout << "Looking for: " << bypass_path << std::endl;

    // Check if file exists
    std::ifstream test_file(bypass_path);
    if (test_file.good()) {
      test_file.close();
      std::cout << "FILE EXISTS - attempting to load..." << std::endl;
      load_vectors_from_csv(bypass_path);
      vectors_loaded = use_preloaded_vectors;
      if (vectors_loaded) {
        std::cout << "SUCCESS: Loaded " << preloaded_vectors.size() << " random vectors from R bypass file" << std::endl;
      } else {
        std::cout << "FAILED: Could not parse random vectors from file" << std::endl;
      }
    } else {
      std::cout << "FILE NOT FOUND - will use C++ random generation" << std::endl;
      std::cout << "To enable bypass: Run R version first to generate this file" << std::endl;
    }
    std::cout << "==================================" << std::endl;
  }
  // Reset vector index on each GetTrace call
  preloaded_vector_idx = 0;

  // Set R's RNG seed (default 10, matches R's GetTrace; overridable via fit.trace_seed)
  if (!use_preloaded_vectors) {
    int seed = getTraceSeedOr(10);
    std::cout << "[GetTrace] trace RNG seed = " << seed
              << (g_trace_seed >= 0 ? " (config override)" : " (builtin default)") << std::endl;
    set_seed(seed);
    GetRNGstate();
  }

  // --- Dimensions & sanity ---
  const int n = Sigma_iX.n_rows;
  const int p = Sigma_iX.n_cols;
  std::cout << "GetTrace: n=" << n << " p=" << p << " nrun=" << nrun << std::endl << std::flush;
  if (Xmat.n_rows != n)          throw std::runtime_error("GetTrace: Xmat.n_rows != Sigma_iX.n_rows");
  if (Xmat.n_cols != p)          throw std::runtime_error("GetTrace: Xmat.n_cols != Sigma_iX.n_cols");
  if ((int)wVec.n_rows != n)     throw std::runtime_error("GetTrace: wVec length != n");
  if (cov1.n_rows != (arma::uword)p || cov1.n_cols != (arma::uword)p)
                                 throw std::runtime_error("GetTrace: cov1 must be p×p");
  if (!Sigma_iX.is_finite() || !Xmat.is_finite() || !wVec.is_finite() ||
      !tauVec.is_finite() || !cov1.is_finite()) {
    throw std::runtime_error("GetTrace: non-finite entries in inputs");
  }

  arma::fmat Sigma_iXt = Sigma_iX.t();   // p×n

  int nrunStart = 0;
  int nrunEnd   = std::max(1, nrun);
  float traceCV = traceCVcutoff + 0.1f;

  arma::fvec tempVec(nrunEnd, arma::fill::zeros);

  while (traceCV > traceCVcutoff) {
    // ensure tempVec has capacity for [nrunStart, nrunEnd)
    if ((int)tempVec.n_rows < nrunEnd) {
      const int old = tempVec.n_rows;
      tempVec.resize(nrunEnd);
      tempVec.rows(old, nrunEnd - 1).zeros();  // zero new tail
    }

    // ---- Phase-2 batched wave: solve all probes of [nrunStart, nrunEnd)
    // with one block-PCG + one batched ψ·U. The probe vectors are generated
    // up front IN THE ORIGINAL ORDER, so the RNG stream (R's rbinom via
    // rademacher_vec) is consumed identically to the serial loop below.
    // SAIGE_NO_BLOCKPCG=1 restores the serial loop.
    if (!isBlockPCGdisabled()) {
      const int nb_cols = nrunEnd - nrunStart;
      arma::fmat Umat(n, nb_cols);
      for (int i = 0; i < nb_cols; ++i)
        Umat.col(i) = rademacher_vec(n);

      arma::fmat Sigma_iU = getPCGofSigmaAndMatrix(wVec, tauVec, Umat,
                                                   maxiterPCG, tolPCG);
      arma::fmat PU = Sigma_iU - Sigma_iX * (cov1 * (Sigma_iXt * Umat));
      arma::fmat AU = getCrossprodMatAndKinMat(Umat);
      if (!AU.is_finite()) throw std::runtime_error("GetTrace: Au non-finite");
      if (!PU.is_finite()) throw std::runtime_error("GetTrace: Pu non-finite");
      for (int i = 0; i < nb_cols; ++i)
        tempVec(nrunStart + i) = arma::dot(AU.col(i), PU.col(i));
    } else
    for (int i = nrunStart; i < nrunEnd; ++i) {
      // uVec: length n, values in {-1, +1}
      if (i == 0 && nrunStart == 0) std::cout << "DEBUG GetTrace: generating rademacher_vec..." << std::endl << std::flush;
      arma::fvec uVec = rademacher_vec(n);
      if (i == 0 && nrunStart == 0) std::cout << "DEBUG GetTrace: calling getPCG1ofSigmaAndVector..." << std::endl << std::flush;

      arma::fvec Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG); // n
      if (i == 0 && nrunStart == 0) std::cout << "DEBUG GetTrace: PCG done, computing Pu..." << std::endl << std::flush;
      arma::fmat cov_safe = cov1;
      arma::fvec Pu = Sigma_iu - Sigma_iX * (cov1 * (Sigma_iX.t() * uVec));
      if (i == 0 && nrunStart == 0) std::cout << "DEBUG GetTrace: calling getCrossprodMatAndKin..." << std::endl << std::flush;
      arma::fvec Au       = getCrossprodMatAndKin(uVec);                                     // n
      if (i == 0 && nrunStart == 0) std::cout << "DEBUG GetTrace: getCrossprodMatAndKin done" << std::endl << std::flush;

#ifdef SAIGE_DEBUG_IO
      // Save C++ Au to file for comparison
      {
        std::string cpp_out = saige_env_path("SAIGE_DEBUG_DIR", "cpp_Au_vec_" + std::to_string(i) + ".txt");
        std::ofstream outfile(cpp_out);
        if (outfile.is_open()) {
          for (int j = 0; j < n; j++) {
            outfile << Au(j) << "\n";
          }
          outfile.close();
        }
      }

      // Compare with R's Au values
      static bool use_r_au = false;  // Set to false to compare, true to override
      {
        std::string filename = "/tmp/r_Au_vec_" + std::to_string(i) + ".txt";
        std::ifstream infile(filename);
        if (infile.is_open()) {
          arma::fvec r_Au(n);
          for (int j = 0; j < n; j++) {
            infile >> r_Au(j);
          }
          infile.close();

          // Compute difference statistics
          arma::fvec diff = Au - r_Au;
          float max_diff = arma::max(arma::abs(diff));
          float mean_diff = arma::mean(diff);
          float norm_diff = arma::norm(diff);
          float norm_Au = arma::norm(Au);
          float norm_r_Au = arma::norm(r_Au);

          if (i < 5) {  // Print for first 5 vectors
            std::cout << "\n=== GRM COMPARISON (vector " << i << ") ===" << std::endl;
            std::cout << "C++ Au[0:5]: " << Au(0) << " " << Au(1) << " " << Au(2) << " " << Au(3) << " " << Au(4) << std::endl;
            std::cout << "R   Au[0:5]: " << r_Au(0) << " " << r_Au(1) << " " << r_Au(2) << " " << r_Au(3) << " " << r_Au(4) << std::endl;
            std::cout << "Max abs diff: " << max_diff << std::endl;
            std::cout << "Mean diff: " << mean_diff << std::endl;
            std::cout << "Norm diff: " << norm_diff << " (relative: " << (norm_diff / norm_r_Au * 100) << "%)" << std::endl;
            std::cout << "|C++ Au|: " << norm_Au << ", |R Au|: " << norm_r_Au << std::endl;
          }

          if (use_r_au) {
            Au = r_Au;  // Use R's Au instead
          }
        }
      }
#endif  // SAIGE_DEBUG_IO
      if ((int)Au.n_rows != n) {
        throw std::runtime_error("GetTrace: Au/Pu bad size");
      }

      if (!Au.is_finite()) {
        throw std::runtime_error("GetTrace: Au non-finite");
      }

      if (!Pu.is_finite()) {
        throw std::runtime_error("GetTrace: Pu non-finite");
      }

      tempVec(i) = arma::dot(Au, Pu);

      // DEBUG: Check magnitude of Au and Pu on first iteration
      if (i == 0 && nrunStart == 0) {
        std::cout << "\n===== C++ GetTrace DEBUG (first vector) =====" << std::endl;
        std::cout << "|Au| (norm): " << arma::norm(Au) << std::endl;
        std::cout << "|Pu| (norm): " << arma::norm(Pu) << std::endl;
        std::cout << "Au[0:5]: " << Au(0) << " " << Au(1) << " " << Au(2) << " " << Au(3) << " " << Au(4) << std::endl;
        std::cout << "Pu[0:5]: " << Pu(0) << " " << Pu(1) << " " << Pu(2) << " " << Pu(3) << " " << Pu(4) << std::endl;
        std::cout << "uVec[0:5]: " << uVec(0) << " " << uVec(1) << " " << uVec(2) << " " << uVec(3) << " " << uVec(4) << std::endl;
        std::cout << "dot(Au, Pu): " << tempVec(i) << std::endl;
        std::cout << "tau[1]: " << tauVec(1) << std::endl;
        std::cout << "=============================================" << std::endl;
      }

      // release temporaries explicitly (optional)
      Au.reset(); Pu.reset(); Sigma_iu.reset(); uVec.reset();
    }

    // Compute CV only on the filled prefix [0, nrunEnd)
    // NOTE: SAIGE uses CV = (sd/mean)/n, not standard CV = sd/mean
    // This is equivalent to the standard error of the mean divided by the mean
    const arma::fvec slice = tempVec.rows(0, nrunEnd - 1);
    const double mu = arma::mean(slice);
    const double sd = arma::stddev(slice);
    // calCV formula: (sd/mean)/veclen
    traceCV = (mu != 0.0) ? static_cast<float>((sd / std::abs(mu)) / nrunEnd)
                          : std::numeric_limits<float>::infinity();

    if (traceCV > traceCVcutoff) {
      std::cerr << "CV for trace random estimator using " << nrunEnd
                << " runs is " << traceCV << " > " << traceCVcutoff
                << " → try " << (nrunEnd + 10) << " runs\n";
      nrunStart = nrunEnd;
      nrunEnd  += 10;
    }
  }

  if (!use_preloaded_vectors) {
    PutRNGstate();
  }

  return arma::mean(tempVec.rows(0, nrunEnd - 1));
}
}

// // INTERNAL: Compute trace of matrix using Monte Carlo estimation
// float GetTrace(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
//   set_seed(200);
//   int Nnomissing = geno.getNnomissing();
//   arma::fmat Sigma_iXt = Sigma_iX.t();
//   arma::fvec Sigma_iu;  
//   arma::fcolvec Pu;
//   arma::fvec Au;
//   arma::fvec uVec;

//   int nrunStart = 0;
//   int nrunEnd = nrun;
//   float traceCV = traceCVcutoff + 0.1;
//   arma::fvec tempVec(nrun);
//   tempVec.zeros();

//   while(traceCV > traceCVcutoff){     
//     //arma::fvec tempVec(nrun);
//     //tempVec.zeros();
//     //arma::fvec tempVec(nrun);
//     //tempVec.zeros();
//     for(int i = nrunStart; i < nrunEnd; i++){
//       Rcpp::NumericVector uVec0;
//       uVec0 = nb(Nnomissing);
//       uVec = as<arma::fvec>(uVec0);
//       uVec = uVec*2 - 1;
//       Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);
//       Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
//       Au = getCrossprodMatAndKin(uVec);
//       tempVec(i) = dot(Au, Pu);
//       Au.clear();
//       Pu.clear();
//       Sigma_iu.clear();
//       uVec.clear();
//     }
//     traceCV = calCV(tempVec);
//     if(traceCV > traceCVcutoff){
//       nrunStart = nrunEnd;
//       nrunEnd = nrunEnd + 10;
//       tempVec.resize(nrunEnd); 
//       cout << "CV for trace random estimator using "<< nrun << " runs is " << traceCV <<  " > " << traceCVcutoff << endl;
//       cout << "try " << nrunEnd << " runs" << endl;      
//     }
//   }

//   float tra = arma::mean(tempVec);
//   tempVec.clear();
//   return(tra);
// }



// Added by SLEE, 04/16/2017
//      This function calculate fixed and random effect coefficients

// R CONNECTION: Core function called from Get_Coef() in R/SAIGE_fitGLMM_fast.R
// Solves for GLMM coefficients using PCG algorithm, returns results to R functions
// Rcpp::List getCoefficients(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG){

//   	int Nnomissing = geno.getNnomissing();
//   	arma::fvec Sigma_iY;
//   	Sigma_iY = getPCG1ofSigmaAndVector(wVec, tauVec, Yvec, maxiterPCG, tolPCG);
//   	int colNumX = Xmat.n_cols;
//   	arma::fmat Sigma_iX(Nnomissing,colNumX);
//   	arma::fvec XmatVecTemp;
//   	for(int i = 0; i < colNumX; i++){
//     		XmatVecTemp = Xmat.col(i);

//     		Sigma_iX.col(i) = getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

//   	}

//   	arma::fmat Xmatt = Xmat.t();
//   	//arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
// 	arma::fmat cov;
// 	try {
// 	  cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
// 	} catch (const std::exception& e) {
// 	  cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
// 	  cout << "inv_sympd failed, inverted with pinv" << endl;
// 	}


//  	arma::fmat Sigma_iXt = Sigma_iX.t();
//   	arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
//   	arma::fvec alpha = cov * SigmaiXtY;

//   	arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
//   	return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
// }



// R CONNECTION: LOCO version called from Get_Coef_LOCO() in R/SAIGE_fitGLMM_fast.R
// Used for leave-one-chromosome-out analysis to avoid genomic inflation
// Rcpp::List getCoefficients_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG){

//         int Nnomissing = geno.getNnomissing();
//         arma::fvec Sigma_iY;

//         Sigma_iY = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Yvec, maxiterPCG, tolPCG);
//         int colNumX = Xmat.n_cols;
//         arma::fmat Sigma_iX(Nnomissing,colNumX);
//         arma::fvec XmatVecTemp;
//         for(int i = 0; i < colNumX; i++){
//                 XmatVecTemp = Xmat.col(i);
//                 Sigma_iX.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

//         }
//         arma::fmat Xmatt = Xmat.t();
//         //arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
//         arma::fmat cov;
//         try {
//           cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
//         } catch (const std::exception& e) {
//           cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
//           cout << "inv_sympd failed, inverted with pinv" << endl;
//         }

//         arma::fmat Sigma_iXt = Sigma_iX.t();
//         arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
//         arma::fvec alpha = cov * SigmaiXtY;

//         arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
//         return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
// }


// 
// Rcpp::List getCoefficients_q_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG){

//         int Nnomissing = geno.getNnomissing();
//         arma::fvec Sigma_iY;

//         Sigma_iY = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Yvec, maxiterPCG, tolPCG);
//         int colNumX = Xmat.n_cols;
//         arma::fmat Sigma_iX(Nnomissing,colNumX);
//         arma::fvec XmatVecTemp;
//         for(int i = 0; i < colNumX; i++){
//                 XmatVecTemp = Xmat.col(i);

//                 Sigma_iX.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

//         }

//         arma::fmat Xmatt = Xmat.t();
//         //arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
//         arma::fmat cov;
//         try {
//           cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
//         } catch (const std::exception& e) {
//           cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
//           cout << "inv_sympd failed, inverted with pinv" << endl;
//         }
//         arma::fmat Sigma_iXt = Sigma_iX.t();
//         arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
//         arma::fvec alpha = cov * SigmaiXtY;

//         arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
//         return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
// }






// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//      This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace

// R CONNECTION: Called from fitglmmaiRPCG() in R for Average Information scoring
// Computes AI matrix and score for variance component estimation
// Rcpp::List getAIScore(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
// arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
// int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){

// 	arma::fmat Sigma_iXt = Sigma_iX.t();
//   	arma::fvec PY1 = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));
// 	//PY1.print("PY1");
//   	arma::fvec APY = getCrossprodMatAndKin(PY1);
// 	//APY.print("APY");
//   	float YPAPY = dot(PY1, APY);
	
//   	float Trace = GetTrace(Sigma_iX, Xmat, wVec, tauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);
// 	//std::cout << "Trace " << Trace << std::endl;
//   	arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG);
// 	//wVec.print("wVec");
// 	//tauVec.print("tauVec");
// 	//PAPY_1.print("PAPY_1");

//   	arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
//   	float AI = dot(APY, PAPY);

//   	return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY1, Named("AI") = AI);
// }

//Rcpp::List fitglmmaiRPCG(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore

// Rcpp::List fitglmmaiRPCG(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec &tauVec,
// arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
// int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){

// 	//double mem1,mem2;
//         //process_mem_usage(mem1, mem2);
//         //std::cout << "fitglmmaiRPCG starts" << std::endl;

// 	//std::cout << "fitglmmaiRPCG starts 2" << std::endl;
//   	Rcpp::List re = getAIScore(Yvec, Xmat,wVec,  tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);

//         //process_mem_usage(mem1, mem2);
//   	float YPAPY = re["YPAPY"];
//   	float Trace = re["Trace"];
//   	float score1 = YPAPY - Trace;
//   	float AI1 = re["AI"];
//   	float Dtau = score1/AI1;
//   	arma::fvec tau0 = tauVec;
//   	tauVec(1) = tau0(1) + Dtau;
//         //std::cout << "fitglmmaiRPCG ends" << std::endl;
//         //std::cout << "AI1: " << AI1 << std::endl;
//         //std::cout << "score1: " << score1 << std::endl;
//         //std::cout << "Dtau: " << Dtau << std::endl;

//   	for(int i=0; i<tauVec.n_elem; ++i) {
//     		if (tauVec(i) < tol){
//       			tauVec(i) = 0;
//     		}
//   	}

//   	float step = 1.0;
//   	while (tauVec(1) < 0.0){

//     		step = step*0.5;
//     		tauVec(1) = tau0(1) + step * Dtau;

//   	}

//   	for(int i=0; i<tauVec.n_elem; ++i) {
//     		if (tauVec(i) < tol){
//       			tauVec(i) = 0;
//     		}
//   	}
//         //process_mem_usage(mem1, mem2);
//    	//std::cout << "VM 3: " << mem1 << "; RSS 3: " << mem2<< std::endl;
//   	return List::create(Named("tau") = tauVec);
// }



/*add for SPA by Wei 04222017*/

// R CONNECTION: Computes Sigma^(-1)*X matrix operations to R functions
// Essential for mixed model fixed effects estimation and coefficient covariance
arma::fmat getSigma_X(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG){


  	int Nnomissing = Xmat.n_rows;
  	int colNumX = Xmat.n_cols;

  	//cout << colNumX << endl;
  	//cout << size(wVec) << endl;
  	//cout << size(tauVec) << endl;


  	arma::fmat Sigma_iX1(Nnomissing,colNumX);
  	arma::fvec XmatVecTemp;

  	for(int i = 0; i < colNumX; i++){
    		XmatVecTemp = Xmat.col(i);
    		Sigma_iX1.col(i) = getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);
  	}
  	return(Sigma_iX1);
}



// R CONNECTION: LOCO version of Sigma^(-1)*X operations to R functions
// Leave-one-chromosome-out computation for unbiased fixed effects estimation
arma::fmat getSigma_X_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;

        //cout << colNumX << endl;
        //cout << size(wVec) << endl;
        //cout << size(tauVec) << endl;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);
        }
        return(Sigma_iX1);
}




// R CONNECTION: Survival-specific Sigma^(-1)*X operations to R functions
// Matrix operations for survival mixed model fixed effects and inference
arma::fmat getSigma_X_Surv(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv,  arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;

        //cout << colNumX << endl;
        //cout << size(wVec) << endl;
        //cout << size(tauVec) << endl;
        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
        arma::fvec x0Vec;


        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                //x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                x0Vec = wVec % XmatVecTemp - sqrtDRN.t() * sqrtDRN * XmatVecTemp;
                //x0Vec.print();

                if(tauVec(1) != 0){

                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv(wVec, tauVec, XmatVecTemp, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
                }else{
                        Sigma_iX1.col(i) = x0Vec;
                }

        }
        return(Sigma_iX1);
}


// R CONNECTION: LOCO survival Sigma^(-1)*X operations to R functions
// Leave-one-chromosome-out matrix operations for survival mixed models
arma::fmat getSigma_X_Surv_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv,  arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;

        //cout << colNumX << endl;
        //cout << size(wVec) << endl;
        //cout << size(tauVec) << endl;
        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
        arma::fvec x0Vec;

        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                //x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                x0Vec = wVec % XmatVecTemp - sqrtDRN.t() * sqrtDRN * XmatVecTemp;

                if(tauVec(1) != 0){

                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_LOCO(wVec, tauVec, XmatVecTemp, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
                }else{
                        Sigma_iX1.col(i) = x0Vec;
                }

        }


        return(Sigma_iX1);
}


// R CONNECTION: Enhanced survival Sigma^(-1)*X operations to R functions
// Optimized matrix operations for survival mixed model coefficient estimation
arma::fmat getSigma_X_Surv_new(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

        //arma::fvec XmatVecTemp1 = Xmat.col(0);
        //arma::fvec Sigma_iX1temp = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, XmatVecTemp1, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, maxiterPCG, tolPCG);
        //cout << "XmatVecTemp1, 1st column in Xmat" << endl;

        arma::fvec x0Vec;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        for(int i = 0; i < colNumX; i++){
                //if (i == 0){
                //      cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
                //}
                XmatVecTemp = Xmat.col(i);
        //      cout << "i is " << i << endl;
        //      cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
                x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);

                if(tauVec(1) != 0){
                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, XmatVecTemp, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
                }else{
                        Sigma_iX1.col(i) = x0Vec;
        //              Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                }
        //      cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;

        }
        return(Sigma_iX1);
}


// R CONNECTION: Enhanced LOCO survival Sigma^(-1)*X operations to R functions
// Optimized leave-one-chromosome-out matrix operations for survival models
arma::fmat getSigma_X_Surv_new_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
          arma::fvec x0Vec;
        //arma::fvec XmatVecTemp1 = Xmat.col(0);
        //arma::fvec Sigma_iX1temp = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, XmatVecTemp1, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG);
        //cout << "XmatVecTemp1, 1st column in Xmat" << endl;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);

        for(int i = 0; i < colNumX; i++){
                //if (i == 0){
                //      cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
                //}
                XmatVecTemp = Xmat.col(i);
                x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
        //      cout << "i is " << i << endl;
        //      cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
                if(tauVec(1)!=0){
                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new_LOCO(wVec, tauVec, XmatVecTemp, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
                }else{
                        //Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                        Sigma_iX1.col(i) = x0Vec;
                }
        //      cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;

        }
        return(Sigma_iX1);
}


arma::fmat  getSigma_X_Surv_new2(arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& Xmat,  arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

        arma::fvec x0Vec;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        for(int i = 0; i < colNumX; i++){
                if (i == 0){
                        cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
                }
                XmatVecTemp = Xmat.col(i);
        //      cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
                x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                cout << "i is " << i << endl;

                if(tauVec(1) != 0){
                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new2(wVec, tauVec, XmatVecTemp, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
                }else{
                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new2(wVec, tauVec, XmatVecTemp, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
                        //Sigma_iX1.col(i) = x0Vec;
        //              Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                }
        //      cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;

        }
        return(Sigma_iX1);
}



arma::fmat  getSigma_X_Surv_new2_LOCO(arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& Xmat,  arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){

        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
          arma::fvec x0Vec;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);

        for(int i = 0; i < colNumX; i++){
                //if (i == 0){
                //      cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
                //}
                XmatVecTemp = Xmat.col(i);
                x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
        //      cout << "i is " << i << endl;
        //      cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
                if(tauVec(1)!=0){
                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_LOCO_new2(wVec, tauVec, XmatVecTemp, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);

                }else{
                        //Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                        Sigma_iX1.col(i) = x0Vec;
                }
        //      cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;

        }
        return(Sigma_iX1);
}






// R CONNECTION: Computes Sigma^(-1)*G vector operations to R functions
// Essential for mixed model genotype effect estimation and association testing
arma::fvec  getSigma_G(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG){
  	arma::fvec Sigma_iG;
  	Sigma_iG = getPCG1ofSigmaAndVector(wVec, tauVec, Gvec, maxiterPCG, tolPCG);
  	return(Sigma_iG);
}


// R CONNECTION: LOCO version of Sigma^(-1)*G operations to R functions
// Leave-one-chromosome-out computation for unbiased genotype effect estimation
arma::fvec  getSigma_G_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Gvec, maxiterPCG, tolPCG);
        return(Sigma_iG);
}




// R CONNECTION: Survival-specific Sigma^(-1)*G operations to R functions
// Vector operations for survival mixed model genotype association testing
arma::fvec  getSigma_G_Surv(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec,  arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;

        arma::fvec x0Vec = wVec % Gvec - sqrtDRN.t() * sqrtDRN * Gvec;

        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv(wVec, tauVec, Gvec, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
        }else{
                //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }

        return(Sigma_iG);
}


// R CONNECTION: LOCO survival Sigma^(-1)*G operations to R functions
// Leave-one-chromosome-out vector operations for survival association testing
arma::fvec  getSigma_G_Surv_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec,  arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;

        arma::fvec x0Vec = wVec % Gvec - sqrtDRN.t() * sqrtDRN * Gvec;

        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_LOCO(wVec, tauVec, Gvec, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
        }else{
                //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }


        return(Sigma_iG);
}


// R CONNECTION: Enhanced survival Sigma^(-1)*G operations to R functions
// Optimized vector operations for survival mixed model genotype testing
arma::fvec  getSigma_G_Surv_new(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        arma::fvec x0Vec;
        x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);


        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG);
        }else{
                //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}



// R CONNECTION: Second-generation survival Sigma^(-1)*G operations to R functions
// Advanced vector operations with enhanced stability for survival association
arma::fvec  getSigma_G_Surv_new2(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
        arma::fvec Sigma_iG;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        arma::fvec x0Vec;
        x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);


        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_new2(wVec, tauVec, Gvec, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
        }else{
        //      //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}


// R CONNECTION: LOCO second-generation survival Sigma^(-1)*G operations to R functions
// Advanced LOCO vector operations for survival mixed model association testing
arma::fvec  getSigma_G_Surv_new2_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
        arma::fvec Sigma_iG;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        arma::fvec x0Vec;
        x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);


        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_LOCO_new2(wVec, tauVec, Gvec, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
        }else{
        //      //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}



// R CONNECTION: Enhanced LOCO survival Sigma^(-1)*G operations to R functions
// Optimized leave-one-chromosome-out vector operations for survival models
arma::fvec  getSigma_G_Surv_new_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
        arma::fvec x0Vec;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);


        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_new_LOCO(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec,diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG);
        }else{
                Sigma_iG = x0Vec;
        //        Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
        }
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}



//This function needs the function getPCG1ofSigmaAndVector and function getCrossprodMatAndKin

arma::fvec GetTrace_q(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
  std::cout << "=== Entering GetTrace_q ===" << std::endl << std::flush;

  // Load precomputed vectors from R if file exists (quantitative trait version)
  static bool load_attempted_q = false;
  if (!load_attempted_q) {
    load_attempted_q = true;
    std::string bypass_path = saige_env_path("SAIGE_BYPASS_DIR", "random_vectors_seed200.csv");
    std::cout << "=== RANDOM VECTOR BYPASS CHECK (GetTrace_q) ===" << std::endl;
    std::cout << "Looking for: " << bypass_path << std::endl;

    // Check if file exists
    std::ifstream test_file(bypass_path);
    if (test_file.good()) {
      test_file.close();
      std::cout << "FILE EXISTS - attempting to load..." << std::endl;
      load_vectors_from_csv(bypass_path);
      if (use_preloaded_vectors) {
        std::cout << "SUCCESS: Loaded " << preloaded_vectors.size() << " random vectors from R bypass file (seed200)" << std::endl;
      } else {
        std::cout << "FAILED: Could not parse random vectors from file" << std::endl;
      }
    } else {
      std::cout << "FILE NOT FOUND - will use C++ random generation" << std::endl;
      std::cout << "To enable bypass: Run R version first to generate this file" << std::endl;
    }
    std::cout << "================================================" << std::endl;
  }
  // Reset vector index on each GetTrace_q call
  preloaded_vector_idx = 0;

  // Set R's RNG seed (default 200, matches R's GetTrace_q; overridable via fit.trace_seed)
  // Then load RNG state so Rf_rbinom() draws from the seeded stream
  if (!use_preloaded_vectors) {
    int seed = getTraceSeedOr(200);
    std::cout << "[GetTrace_q] trace RNG seed = " << seed
              << (g_trace_seed >= 0 ? " (config override)" : " (builtin default)") << std::endl;
    set_seed(seed);
    GetRNGstate();
  }

  const int n = Sigma_iX.n_rows;
  arma::fmat Sigma_iXt = Sigma_iX.t();

  int nrunStart = 0;
  int nrunEnd = std::max(1, nrun);
  float traceCV  = traceCVcutoff + 0.1f;
  float traceCV0 = traceCVcutoff + 0.1f;

  arma::fvec tempVec(nrunEnd, arma::fill::zeros);
  arma::fvec tempVec0(nrunEnd, arma::fill::zeros);

  while ((traceCV > traceCVcutoff) || (traceCV0 > traceCVcutoff)) {
    if ((int)tempVec.n_rows < nrunEnd) {
      const int old = tempVec.n_rows;
      tempVec.resize(nrunEnd);
      tempVec.rows(old, nrunEnd - 1).zeros();
      tempVec0.resize(nrunEnd);
      tempVec0.rows(old, nrunEnd - 1).zeros();
    }

    // ---- Phase-2 batched wave (see GetTrace): one block-PCG + one batched
    // ψ·U per wave; RNG stream consumed in the original per-probe order.
    if (!isBlockPCGdisabled()) {
      const int nb_cols = nrunEnd - nrunStart;
      arma::fmat Umat(n, nb_cols);
      for (int i = 0; i < nb_cols; ++i)
        Umat.col(i) = rademacher_vec(n);

      arma::fmat Sigma_iU = getPCGofSigmaAndMatrix(wVec, tauVec, Umat,
                                                   maxiterPCG, tolPCG);
      arma::fmat PU = Sigma_iU - Sigma_iX * (cov1 * (Sigma_iXt * Umat));
      arma::fmat AU = getCrossprodMatAndKinMat(Umat);
      for (int i = 0; i < nb_cols; ++i) {
        tempVec(nrunStart + i)  = arma::dot(AU.col(i), PU.col(i));
        tempVec0(nrunStart + i) = arma::dot(Umat.col(i), PU.col(i));
      }
    } else
    for (int i = nrunStart; i < nrunEnd; ++i) {
      arma::fvec uVec = rademacher_vec(n);

      arma::fvec Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);
      arma::fvec Pu = Sigma_iu - Sigma_iX * (cov1 * (Sigma_iXt * uVec));
      arma::fvec Au = getCrossprodMatAndKin(uVec);

      tempVec(i)  = arma::dot(Au, Pu);   // trace for kinship component (tau[1])
      tempVec0(i) = arma::dot(uVec, Pu); // trace for identity component (tau[0])

      Au.reset(); Pu.reset(); Sigma_iu.reset(); uVec.reset();
    }

    // Compute CV for both trace estimators
    {
      const arma::fvec slice = tempVec.rows(0, nrunEnd - 1);
      const double mu = arma::mean(slice);
      const double sd = arma::stddev(slice);
      traceCV = (mu != 0.0) ? static_cast<float>((sd / std::abs(mu)) / nrunEnd)
                             : std::numeric_limits<float>::infinity();
    }
    {
      const arma::fvec slice0 = tempVec0.rows(0, nrunEnd - 1);
      const double mu0 = arma::mean(slice0);
      const double sd0 = arma::stddev(slice0);
      traceCV0 = (mu0 != 0.0) ? static_cast<float>((sd0 / std::abs(mu0)) / nrunEnd)
                               : std::numeric_limits<float>::infinity();
    }

    if ((traceCV > traceCVcutoff) || (traceCV0 > traceCVcutoff)) {
      std::cout << "CV for trace random estimator using " << nrunEnd
                << " runs is " << traceCV << " / " << traceCV0
                << " (> " << traceCVcutoff << ")" << std::endl;
      nrunStart = nrunEnd;
      nrunEnd += 10;
      std::cout << "try " << nrunEnd << " runs" << std::endl;
    }
  }

  arma::fvec traVec(2);
  traVec(1) = arma::mean(tempVec.rows(0, nrunEnd - 1));
  traVec(0) = arma::mean(tempVec0.rows(0, nrunEnd - 1));
  // Restore R's RNG state after all vectors generated
  if (!use_preloaded_vectors) {
    PutRNGstate();
  }

  std::cout << "GetTrace_q: Trace[0] (identity) = " << traVec(0)
            << ", Trace[1] (kinship) = " << traVec(1) << std::endl;
  return traVec;
}

//Rcpp::List getAIScore_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){


//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace

// Rcpp::List getAIScore_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
// arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
// int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){


//   	arma::fmat Sigma_iXt = Sigma_iX.t();
//   	arma::fmat Xmatt = Xmat.t();

//   	//arma::fmat cov1 = inv_sympd(Xmatt * Sigma_iX);
//         arma::fmat cov1;
//         try {
//           cov1 = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
//         } catch (const std::exception& e) {
//           cov1 = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
//           cout << "inv_sympd failed, inverted with pinv" << endl;
//         }


//   	arma::fvec PY1 = Sigma_iY - Sigma_iX * (cov1 * (Sigma_iXt * Yvec));
//   	arma::fvec APY = getCrossprodMatAndKin(PY1);

//   	float YPAPY = dot(PY1, APY);

//   	arma::fvec A0PY = PY1; ////Quantitative


//   	float YPA0PY = dot(PY1, A0PY); ////Quantitative


//   	arma::fvec Trace = GetTrace_q(Sigma_iX, Xmat, wVec, tauVec, cov1, nrun, maxiterPCG, tolPCG, traceCVcutoff);

//   	arma::fmat AI(2,2);
//   	arma::fvec PA0PY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, A0PY, maxiterPCG, tolPCG);
//   	arma::fvec PA0PY = PA0PY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PA0PY_1));

//   	AI(0,0) =  dot(A0PY, PA0PY);

//   	//cout << "A1(0,0) " << AI(0,0)  << endl;
//   	arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG);
//   	arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PAPY_1));
//   	AI(1,1) = dot(APY, PAPY);

//   	AI(0,1) = dot(A0PY, PAPY);

//   	AI(1,0) = AI(0,1);

//   	//cout << "AI " << AI << endl;
//   	//cout << "Trace " << Trace << endl;
//   	//cout << "YPAPY " << YPAPY << endl;
//   	//cout << "cov " << cov1 << endl;
// 	return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("YPA0PY") = YPA0PY,Named("Trace") = Trace,Named("PY") = PY1,Named("AI") = AI);

// }






//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace

// Rcpp::List getAIScore_q_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){

//         int Nnomissing = geno.getNnomissing();
//         arma::fvec Sigma_iY1;
//         Sigma_iY1 = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Yvec, maxiterPCG, tolPCG);

// //	for(int j = 0; j < 10; j++){
// //                std::cout << "Sigma_iY1(j): " << Sigma_iY1(j) << std::endl;
// //        }


//         int colNumX = Xmat.n_cols;
//         arma::fmat Sigma_iX1(Nnomissing,colNumX);
//         arma::fvec XmatVecTemp;

//         for(int i = 0; i < colNumX; i++){
//                 XmatVecTemp = Xmat.col(i);

//                 Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

//         }


//         //rma::fmat Sigma_iX1t = Sigma_iX1.t();
//         arma::fmat Xmatt = Xmat.t();

//         //arma::fmat cov1 = inv_sympd(Xmatt * Sigma_iX1);
//         arma::fmat cov1;
//         try {
//           cov1 = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX1));
//         } catch (const std::exception& e) {
//           cov1 = arma::pinv(arma::symmatu(Xmatt * Sigma_iX1));
//           cout << "inv_sympd failed, inverted with pinv" << endl;
//         }


//         //cout << "cov " << cov1 << endl;


// 	return Rcpp::List::create(Named("cov") = cov1, Named("Sigma_iX") = Sigma_iX1, Named("Sigma_iY") = Sigma_iY1);
//         //return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace,Named("Sigma_iY") = Sigma_iY1, Named("Sigma_iX") = Sigma_iX1, Named("PY") = PY1, Named("AI") = AI, Named("cov") = cov1);
// }



//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore_q

// Rcpp::List fitglmmaiRPCG_q_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){

//         arma::uvec zeroVec = (tauVec < tol); //for Quantitative, GMMAT
//         Rcpp::List re = getAIScore_q_LOCO(Yvec, Xmat, wVec, tauVec, nrun, maxiterPCG, tolPCG, traceCVcutoff);
// //return Rcpp::List::create(Named("cov") = cov1, Named("Sigma_iX") = Sigma_iX1, Named("Sigma_iY") = Sigma_iY1);
//         arma::fmat cov = re["cov"];
//         arma::fmat Sigma_iX = re["Sigma_iX"];
//         arma::fmat Sigma_iXt = Sigma_iX.t();

//         arma::fvec alpha1 = cov * (Sigma_iXt * Yvec);
//         arma::fvec Sigma_iY = re["Sigma_iY"];
//         arma::fvec eta1 = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha1) / wVec;
// 	return List::create(Named("tau") = tauVec, Named("cov") = cov, Named("alpha") = alpha1, Named("eta") = eta1);
// }



//Rcpp::List fitglmmaiRPCG_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){



//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore_q

// Rcpp::List fitglmmaiRPCG_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec &tauVec,
// arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
// int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){

//   	arma::uvec zeroVec = (tauVec < tol); //for Quantitative, GMMAT
// 	Rcpp::List re = getAIScore_q(Yvec, Xmat,wVec,  tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);

//   	float YPAPY = re["YPAPY"];
//   	float YPA0PY = re["YPA0PY"]; //for Quantitative
//   	arma::fvec Trace = re["Trace"]; //for Quantitative

//   	float score0 = YPA0PY - Trace(0); //for Quantitative
//   	float score1 = YPAPY - Trace(1); //for Quantitative
//   	arma::fvec scoreVec(2); //for Quantitative
//   	scoreVec(0) = score0; //for Quantitative
//   	scoreVec(1) = score1; //for Quantitative

//     //for Quantitative
//   	arma::fmat AI = re["AI"];
//   	//cout << "0,0" << AI(0,0) << endl;
//   	//cout << "0,1" << AI(0,1) << endl;
//   	//cout << "1,0" << AI(1,0) << endl;
//   	//cout << "1,1" << AI(1,1) << endl;

//   	arma::fvec Dtau = solve(AI, scoreVec);


//   	arma::fvec tau0 = tauVec;
//   	tauVec = tau0 + Dtau;


//   	tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); //for Quantitative Copied from GMMAT  

//   	float step = 1.0;


//   	//cout << "tau2 " << tauVec(0) << " " << tauVec(1) << endl;
//   	while (tauVec(0) < 0.0 || tauVec(1)  < 0.0){ //for Quantitative
//      	//	cout << "tauVec Here: " << tauVec << endl;
//     		step = step*0.5;
//     		tauVec = tau0 + step * Dtau; //for Quantitative
//     	//	cout << "tau_4: " << tauVec << endl;
//     		tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); //for Quantitative Copied from GMMAT
//     	//	cout << "tau_5: " << tauVec << endl;
//  	}



//   	tauVec.elem( find(tauVec < tol) ).zeros();
// 	return List::create(Named("tau") = tauVec);

//   	//return List::create(Named("tau") = tauVec, Named("cov") = cov, Named("alpha") = alpha1, Named("eta") = eta1);
// }



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_usingSubMarker : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_M_Submarker;
        unsigned int m_M;
        arma::ivec subMarkerIndex ;

        // product that I have accumulated
        arma::fvec m_bout;


        // constructors
        CorssProd_usingSubMarker(arma::fcolvec & y)
                : m_bVec(y) {

                //m_Msub = geno.getMsub();
                subMarkerIndex = getSubMarkerIndex();
                m_M_Submarker = subMarkerIndex.n_elem;
                m_N = geno.getNnomissing();
                m_bout.zeros(m_N);
        }
        CorssProd_usingSubMarker(const CorssProd_usingSubMarker& CorssProd_usingSubMarker, Split)
                : m_bVec(CorssProd_usingSubMarker.m_bVec)
        {

                m_N = CorssProd_usingSubMarker.m_N;
                //m_M = CorssProd_usingSubMarker.m_M;
                m_M_Submarker = CorssProd_usingSubMarker.m_M_Submarker;
                subMarkerIndex = CorssProd_usingSubMarker.subMarkerIndex;
                m_bout.zeros(m_N);

        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                float val1;
                int j;
                for(unsigned int i = begin; i < end; i++){
                        j = subMarkerIndex[i];
//			std::cout << "j: " << j << std::endl;	
                        geno.Get_OneSNP_StdGeno(j, &vec);
                        val1 = dot(vec,  m_bVec);
                        m_bout += val1 * (vec);
                }
        }

        // join my value with that of another InnerProduct
        void join(const  CorssProd_usingSubMarker & rhs) {
        m_bout += rhs.m_bout;
        }
};



arma::fvec parallelCrossProd_usingSubMarker(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        int m_M_Submarker = getSubMarkerNum();

//	std::cout << "m_M_Submarker: " << m_M_Submarker << std::endl;
        CorssProd_usingSubMarker CorssProd_usingSubMarker(bVec);
//	std::cout << "m_M_Submarker: 2 " << m_M_Submarker << std::endl;
  // call paralleReduce to start the work
        parallelReduce(0, m_M_Submarker, CorssProd_usingSubMarker);
//	std::cout << "m_M_Submarker: 3 " << m_M_Submarker << std::endl;
//	std::cout << "CorssProd_usingSubMarker.m_bout " << CorssProd_usingSubMarker.m_bout << std::endl;
  // return the computed product
        //cout << "Msub: " << Msub << endl;
        //for(int i=0; i<100; ++i)
        //{
        //      cout << (CorssProd_usingSubMarker.m_bout/m_M_Submarker)[i] << ' ';
        //}
//        cout << endl;

//	cout << (CorssProd_usingSubMarker.m_bout).n_elem << endl;
        return CorssProd_usingSubMarker.m_bout/m_M_Submarker;
}




arma::fvec getCrossprodMatAndKin_usingSubMarker(arma::fcolvec& bVec){

        arma::fvec crossProdVec = parallelCrossProd_usingSubMarker(bVec) ;

        return(crossProdVec);
}









//std::vector<int> calGRMvalueUsingSubMarker_forOneInv(int sampleIndex, float relatednessCutoff){
//        //sampleIndex starts with 0
//        std::vector<int> relatedIndex;
//        int Ntotal = geno.getNnomissing();
//        arma::fcolvec bindexvec(Ntotal);
//        bindexvec.zeros();
//        bindexvec(sampleIndex) = 1;
//        arma::fvec crossProdVec = getCrossprodMatAndKin_usingSubMarker(bindexvec);
//        for(int i=sampleIndex; i< Ntotal; i++){
//                if(crossProdVec(i) >= relatednessCutoff){
//                        relatedIndex.push_back(i);
//                }
//        }
//        return(relatedIndex);
//}




//The code below is from http://gallery.rcpp.org/articles/parallel-inner-product/
struct InnerProduct : public Worker
{
   // source vectors
   std::vector<float> x;
   std::vector<float> y;

   // product that I have accumulated
   float product;

   // constructors
   InnerProduct(const std::vector<float> x, const std::vector<float> y)
      : x(x), y(y), product(0) {}
   InnerProduct(const InnerProduct& innerProduct, Split)
      : x(innerProduct.x), y(innerProduct.y), product(0) {}

   // process just the elements of the range I've been asked to
   void operator()(std::size_t begin, std::size_t end) {
      product += std::inner_product(x.begin() + begin,
                                    x.begin() + end,
                                    y.begin() + begin,
                                    0.0);
   }

   // join my value with that of another InnerProduct
   void join(const InnerProduct& rhs) {
     product += rhs.product;
   }
};



// R CONNECTION: Computes parallel inner product of two vectors to R functions
// High-performance dot product calculation using parallel processing for large vectors
float parallelInnerProduct(std::vector<float> &x, std::vector<float> &y) {

   int xsize = x.size();
   // declare the InnerProduct instance that takes a pointer to the vector data
   InnerProduct innerProduct(x, y);

   // call paralleReduce to start the work
   parallelReduce(0, x.size(), innerProduct);

   // return the computed product
   return innerProduct.product/xsize;
}



// R CONNECTION: Calculates GRM value for a specific sample pair to R functions
// Computes genomic relationship between two samples for kinship analysis
float calGRMValueforSamplePair(arma::ivec &sampleidsVec){
        //std::vector<float> stdGenoforSamples = geno.Get_Samples_StdGeno(sampleidsVec);
        geno.Get_Samples_StdGeno(sampleidsVec);
	//std::cout << "here5" << std::endl;
	//for(int i = 0; i < 10; i++){
	//	std::cout << geno.stdGenoforSamples[i] << " ";
	//}
	//std::cout << std::endl;
	//std::cout << geno.stdGenoforSamples.size() << std::endl;
        int Ntotal = geno.getNnomissing();
        float grmValue;
	std::vector<float> stdGenoforSamples2;
	//std::cout << "here5b" << std::endl;
	//std::cout << sampleidsVec.n_elem << std::endl;
	//std::cout << "here5c" << std::endl;
        if(sampleidsVec.n_elem == 2){
                std::vector<float> s1Vec;
                //s1Vec.zeros(Ntotal);

                std::vector<float> s2Vec;
                //arma::fvec s2Vec;
                //s2Vec.zeros(Ntotal);

                for(int i = 0; i < Ntotal; i++){
                        //s1Vec[i] = stdGenoforSamples[i*2+0];
                        s1Vec.push_back(geno.stdGenoforSamples[i*2]);
                        s2Vec.push_back(geno.stdGenoforSamples[i*2+1]);
                }
                grmValue = parallelInnerProduct(s1Vec, s2Vec);
                //grmValue = innerProductFun(s1Vec, s2Vec);
        }else{
	//	std::cout << "here5d" << std::endl;
	//	std::cout << "geno.stdGenoforSamples.size() " << geno.stdGenoforSamples.size() << std::endl;
		stdGenoforSamples2.clear();
		for (int i=0; i< geno.stdGenoforSamples.size(); i++){
			//std::cout << i << " " << geno.stdGenoforSamples[i] << " ";
        		stdGenoforSamples2.push_back(geno.stdGenoforSamples[i]);
		}
	//	std::cout << std::endl;
	//	std::cout << "here6" << std::endl;
                grmValue = parallelInnerProduct(stdGenoforSamples2, geno.stdGenoforSamples);
                //grmValue = innerProductFun(stdGenoforSamples2, geno.stdGenoforSamples);
	//	std::cout << "here7" << std::endl;
        }
        return(grmValue);
}


//Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
//arma::sp_fmat createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){




// R CONNECTION: Creates sparse kinship matrix from marker subset to R functions
// Constructs efficient sparse representation of genomic relationships for mixed models
// Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){

//         int nSubMarker = markerIndexVec.n_elem;
//         int Ntotal = geno.getNnomissing();
//         std::vector<unsigned int>     iIndexVec;
//         std::vector<unsigned int>     iIndexVec2;
//         std::vector<unsigned int>     jIndexVec;
//         std::vector<unsigned int>     jIndexVec2;
//         std::vector<unsigned int>     allIndexVec;
//         std::vector<float>     kinValueVec;
//         std::vector<float>     kinValueVec2;
// 	std::vector<float> stdGenoMultiMarkers;	
// 	stdGenoMultiMarkers.resize(Ntotal*nSubMarker);

// 	//std::cout << "createSparseKin1" << std::endl;
// 	size_t sizeTemp;
// 	float kinValue;
// 	float kinValueTemp;
// 	//std::cout << "createSparseKin1b" << std::endl;

// 	Get_MultiMarkersBySample_StdGeno(markerIndexVec, stdGenoMultiMarkers);
// 	std::cout << "createSparseKin2" << std::endl;
// 	//arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), Ntotal, nSubMarker);
// 	arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), nSubMarker, Ntotal);
// 	//std::cout << "createSparseKin3" << std::endl;
// 	//std::cout << "stdGenoMultiMarkersMat.n_rows: " << stdGenoMultiMarkersMat.n_rows << std::endl;
// 	//std::cout << "stdGenoMultiMarkersMat.n_cols: " << stdGenoMultiMarkersMat.n_cols << std::endl;



//         for(unsigned int i=0; i< Ntotal; i++){
//               for(unsigned int j = i; j < Ntotal; j++){
//                         //kinValueTemp = arma::dot(stdGenoMultiMarkersMat.row(i), stdGenoMultiMarkersMat.row(j));
// 			if(j > i){
//                 		kinValueTemp = arma::dot(stdGenoMultiMarkersMat.col(i), stdGenoMultiMarkersMat.col(j));
//                 		kinValueTemp = kinValueTemp/nSubMarker;
//                 		if(kinValueTemp >= relatednessCutoff){
// //                              if(i == 0){
//                                 //std::cout << "kinValueTemp: " << kinValueTemp << std::endl;
//                                 //std::cout << "relatednessCutoff: " << relatednessCutoff << std::endl;
//                                 //std::cout << "i: " << i << std::endl;
// //                              std::cout << "j: " << j;
// //                              }
//                         		iIndexVec.push_back(i);
// 					jIndexVec.push_back(j);

//                 		}
// 			}else{
// 				iIndexVec.push_back(i);
// 				jIndexVec.push_back(j);
// 			}
//         	}
// 	}
	
// 	arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
//         size_t ni = iIndexVec.size();
//         kinValueVec.resize(ni);
//         std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

//         int Mmarker = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
// 		//geno.getM();
//         for(size_t i=0; i< Mmarker; i++){
//                 geno.Get_OneSNP_StdGeno(i, temp);
//                 for(size_t j=0; j < ni; j++){
//                         kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
//                 }

//         }






// /*
// //	(stdGenoMultiMarkersMat.row(0)).print("stdGenoMultiMarkersMat.row(0):");
// 	//std::cout << stdGenoMultiMarkersMat << std::endl;
// 	//std::cout << stdGenoMultiMarkersMat.row(487) << std::endl;	
// 	omp_set_dynamic(0);     // Explicitly disable dynamic teams
//         omp_set_num_threads(16); // Use 16 threads for all consecutive parallel regions
// 	int totalCombination = Ntotal*(Ntotal-1)/2 - 1;

// 	#pragma omp parallel
// 	{
// 	std::vector<unsigned int> vec_privatei;	
// 	std::vector<unsigned int> vec_privatej;	
// 	#pragma omp for nowait //fill vec_private in parallel
// 	for(int k = 0; k < totalCombination; k++){
//         	int i = k / Ntotal;
//         	int j = k % Ntotal;
//         	if((j <= i)){
//             		i = Ntotal - i - 2;
//             		j = Ntotal - j - 1;
//         	}

// //        for(i=0; i< Ntotal; i++){
// //		for(j = i; j < Ntotal; j++){
// 			//kinValueTemp = arma::dot(stdGenoMultiMarkersMat.row(i), stdGenoMultiMarkersMat.row(j));
// 		kinValueTemp = arma::dot(stdGenoMultiMarkersMat.col(i), stdGenoMultiMarkersMat.col(j));
// 		kinValueTemp = kinValueTemp/nSubMarker;
// 		if(kinValueTemp >= relatednessCutoff){
// //				if(i == 0){
// 				//std::cout << "kinValueTemp: " << kinValueTemp << std::endl;
// 				//std::cout << "relatednessCutoff: " << relatednessCutoff << std::endl;
// 				//std::cout << "i: " << i << std::endl;
// //				std::cout << "j: " << j;
// //				}
// 			vec_privatei.push_back((unsigned int)i);
// 			//allIndexVec.push_back(i);
// 			//iIndexVec.push_back(i);
// 			//iIndexVec.push_back(i);
// 			//allIndexVec.push_back(j);
// 			vec_privatej.push_back((unsigned int)j);
				
								
// 		}
// 	}
// //	#pragma omp critical
// 	#pragma omp for schedule(static) ordered
//     	for(int i=0; i<omp_get_num_threads(); i++) {
//         	#pragma omp ordered
//         	iIndexVec.insert(iIndexVec.end(), vec_privatei.begin(), vec_privatei.end());  
//         	jIndexVec.insert(jIndexVec.end(), vec_privatej.begin(), vec_privatej.end());  
//     	}
// //		}
// 	}
// //	int nall = allIndexVec.size();
// //	 std::cout << "nall: " << nall << std::endl;
// //	int k = 0;
// //	while(k < nall){
// 	//	std::cout << "allIndexVec[k]: " << k << " " << allIndexVec[k] << std::endl;
// 	//	std::cout << "allIndexVec[k+1]: " << k+1 << " " << allIndexVec[k+1] << std::endl;
// //        	iIndexVec.push_back(allIndexVec[k]);
// //                jIndexVec.push_back(allIndexVec[k+1]);
// //		k = k + 2;
// //        }
// //	allIndexVec.clear();

// 	for(int k = 0; k < Ntotal; k++){
// 		iIndexVec.push_back((unsigned int)k);
// 		jIndexVec.push_back((unsigned int)k);
// 	}

//         arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
//         size_t ni = iIndexVec.size();
//         //size_t ni = nall/2 + Ntotal;
//         kinValueVec.resize(ni);
//         std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

//         int Mmarker = geno.getM();
//         for(size_t i=0; i< Mmarker; i++){
//                 geno.Get_OneSNP_StdGeno(i, temp);
//                 for(size_t j=0; j < ni; j++){
// //                for(size_t k=0; k < nall/2; k++){
//                         kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
// //                        kinValueVec[j] = kinValueVec[j] + (((*temp)[allIndexVec[k*2]])*((*temp)[allIndexVec[k*2+1]]))/Mmarker;
//                 }
// //		for(size_t k=nall/2; k < ni; k++){
			
// //			kinValueVec[j] = kinValueVec[j] + (((*temp)[allIndexVec[k*2]])*((*temp)[allIndexVec[k*2+1]]))/Mmarker;

// //		}
//         }	


// */   // end of the openMP version 

// 	std::cout << "ni: " << ni << std::endl;
// /*	for(size_t j=0; j < 10; j++){
// 		std::cout << "iIndexVec[j]: " << iIndexVec[j] << std::endl;
// 		std::cout << "jIndexVec[j]: " << jIndexVec[j] << std::endl;
// 		std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
// 	}
// */
// 	for(size_t j=0; j < ni; j++){
// 		if(kinValueVec[j] >= relatednessCutoff){
// 	//	std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
// 			kinValueVec[j] = tauVec(1)*kinValueVec[j];
// 			iIndexVec2.push_back(iIndexVec[j]+1);
// 			jIndexVec2.push_back(jIndexVec[j]+1);
// 			if(iIndexVec[j] == jIndexVec[j]){
// 				kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));	
// 			}
// 			kinValueVec2.push_back(kinValueVec[j]);
// 		}

// 	}

// //	std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;

// 	//arma::fvec x(kinValueVec2);
// 	//arma::umat locations(iIndexVec2);
// 	//arma::uvec jIndexVec2_b(jIndexVec2);
// 	//locations.insert_cols(locations.n_cols, jIndexVec2_b); 
// 	//arma::umat locationst = locations.t();
// 	//locations.clear();
	
// 	//create a sparse Sigma
// //	arma::sp_fmat sparseSigma(locationst, x);
// //	arma::sp_fmat sparseSigmab  = arma::symmatu(sparseSigma);
// 	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);
// //	return sparseSigmab;
// }




arma::fmat getColfromStdGenoMultiMarkersMat(arma::uvec & a){
	return((geno.stdGenoMultiMarkersMat).cols(a));
}


int getNColStdGenoMultiMarkersMat(){
	return((geno.stdGenoMultiMarkersMat).n_cols);
}


int getNRowStdGenoMultiMarkersMat(){
        return((geno.stdGenoMultiMarkersMat).n_rows);
}



// R CONNECTION: Sets subset marker indices for sparse GRM construction from R functions
// Configures which markers to use for sparse kinship matrix computation
void setSubMarkerIndex(arma::ivec &subMarkerIndexRandom){
	geno.subMarkerIndex = subMarkerIndexRandom;
//	std::cout << "(geno.subMarkerIndex).n_elem: " << (geno.subMarkerIndex).n_elem << std::endl;
	int Nnomissing = geno.getNnomissing();
	(geno.stdGenoMultiMarkersMat).set_size(subMarkerIndexRandom.n_elem, Nnomissing);
}


// R CONNECTION: Sets relatedness threshold for kinship matrix filtering from R functions
// Configures minimum genetic similarity required to retain sample relationships
void setRelatednessCutoff(float a){
	geno.relatednessCutoff = a;
}



// REMOVED: innerProduct() - use getInnerProd() from src/UTIL.cpp instead


//Rcpp::List refineKin(std::vector<unsigned int> &iIndexVec, std::vector<unsigned int> & jIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
//Rcpp::List refineKin(arma::imat &iMat, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){


// R CONNECTION: Refines kinship matrix by applying relatedness threshold to R functions
// Filters and optimizes kinship relationships based on genetic similarity cutoff
// Rcpp::List refineKin(float relatednessCutoff){
//         std::vector<unsigned int>     iIndexVec2;
//         std::vector<unsigned int>     jIndexVec2;
// //	std::vector<float>     kinValueVec;
//         std::vector<float>     kinValueVec2;
//  //       std::vector<float>     kinValueVec_orig; //for test original kinship

// 	arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
// 	(*temp).clear();
//         //size_t ni = iIndexVec.size();
//         //size_t ni = iMat.n_rows;
//         size_t ni = geno.indiceVec.size();
// 	std::cout << "ni: " << ni << std::endl;
 
// 	initKinValueVecFinal(ni);

// //	std::cout << "OKK: "  << std::endl;
// //        kinValueVec.resize(ni);
// //        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

//         //int Mmarker = geno.getM();
//         int Mmarker = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM(); 

//         //for(size_t i=0; i< Mmarker; i++){
//         //        geno.Get_OneSNP_StdGeno(i, temp);
//         //        for(size_t j=0; j < ni; j++){
//         //                kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
//         //        }
//         //}
// 	//arma::fvec kinValueVecTemp;
// 	arma::fvec kinValueVecTemp2;
// 	arma::fvec GRMvec;
// 	GRMvec.set_size(ni);
// 	//int Mmarker_mafgr1perc = 0;
//   	for(size_t i=0; i< Mmarker; i++){
// //		std::cout << "OKKK: "  << std::endl;
// //		std::cout << "Mmarker: " << std::endl;

// //                geno.Get_OneSNP_StdGeno(i, temp);
// 		float freqv = geno.alleleFreqVec[i];
// 		//if(freqv >= minMAFtoConstructGRM && freqv <= 1-minMAFtoConstructGRM){
// 		//Mmarker_mafgr1perc = Mmarker_mafgr1perc + 1;

//                 geno.Get_OneSNP_Geno(i);
// 		float invstdv = geno.invstdvVec[i];
// 		geno.setSparseKinLookUpArr(freqv, invstdv);			

// 		//std::cout << "freqv: " << freqv << std::endl;
// 		//std::cout << "invstdv: " << invstdv << std::endl;
// 		//for (int j = 0; j < 3; j++){
// 		//	std::cout << geno.sKinLookUpArr[j][0] << std::endl;	
// 		//	std::cout << geno.sKinLookUpArr[j][1] << std::endl;	
// 		//	std::cout << geno.sKinLookUpArr[j][2] << std::endl;	

// 		//}
// 		//std::cout << "geno.m_OneSNP_StdGeno(i) " << geno.m_OneSNP_StdGeno(i) <<  std::endl;	
// 		//kinValueVecTemp = parallelcalsparseGRM(iMat);
// //		parallelcalsparseGRM(iMat, GRMvec);

// 		parallelcalsparseGRM(GRMvec);
// 		//std::cout << "kinValueVecTemp.n_elem: " << kinValueVecTemp.n_elem << std::endl;
// //		std::cout << "OKKK2: "  << std::endl;
// 		parallelsumTwoVec(GRMvec);
// //		for(size_t j=0; j< ni; j++){
// //			(geno.kinValueVecFinal)[j] = (geno.kinValueVecFinal)[j] + GRMvec(j);
// //		}
// 		(*temp).clear();
// 	   //}//if(freqv >= 0.01 && freqv <= 0.99){
// 		//kinValueVecTemp.clear();
//         }



//        // for(size_t j=0; j < 100; j++){
//        //         std::cout << "iIndexVec[j]: " << iIndexVec[j] << std::endl;
//        //         std::cout << "jIndexVec[j]: " << jIndexVec[j] << std::endl;
//        //         std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
//        // }

// 	int a1;
// 	int a2;
//         for(size_t j=0; j < ni; j++){
// 		geno.kinValueVecFinal[j] = (geno.kinValueVecFinal[j]) /(Mmarker);

// //		std::cout << "j: " << j << " geno.kinValueVecFinal[j]: " << geno.kinValueVecFinal[j] << std::endl;
//             //    if(geno.kinValueVecFinal[j] >= relatednessCutoff){
//                 if((geno.kinValueVecFinal[j]) >= relatednessCutoff){
//         //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
// 			//kinValueVec_orig.push_back((geno.kinValueVecFinal)[j]); //for test	
//                         //(geno.kinValueVecFinal)[j] = tauVec(1)*(geno.kinValueVecFinal)[j];
//                         //(geno.kinValueVecFinal)[j] = tauVec(1)*(geno.kinValueVecFinal)[j];
//  				 a1 = (geno.indiceVec)[j].first + 1;
// 				 a2 = (geno.indiceVec)[j].second + 1;
// 				 iIndexVec2.push_back(a1);
// 				 jIndexVec2.push_back(a2);

//                         kinValueVec2.push_back((geno.kinValueVecFinal)[j]);
//                 }

//         }



// 	std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
// 	//return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2,  Named("kinValue_orig") = kinValueVec_orig);	
// 	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);	
// }



// R CONNECTION: Shortens kinship list by removing low-relatedness pairs to R functions
// Optimizes kinship matrix storage by filtering out weakly related sample pairs
// Rcpp::List shortenList(arma::imat &iMat, arma::fvec &kinValueVecTemp, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
// 	        std::vector<unsigned int>     iIndexVec2;
//         std::vector<unsigned int>     jIndexVec2;
// 	std::vector<float>     kinValueVec2;
// 	size_t ni = iMat.n_rows;

// 	for(size_t j=0; j < ni; j++){
//                 if(kinValueVecTemp(j) >= relatednessCutoff){
//         //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
//                         kinValueVecTemp(j) = tauVec(1)*(kinValueVecTemp(j));
//                         iIndexVec2.push_back(iMat(j,1)+1);
//                         //iIndexVec2.push_back(iIndexVec[j]+1);
//                         jIndexVec2.push_back(iMat(j,2)+1);
//                         //jIndexVec2.push_back(jIndexVec[j]+1);
//         //                if(iIndexVec[j] == jIndexVec[j]){
//         //                        kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));
//         //                }

//                         if(iMat(j,1) == iMat(j,2)){
//                                 kinValueVecTemp(j) = kinValueVecTemp(j) + tauVec(0)/(wVec(iMat(j,1)));
//                         }

//                         kinValueVec2.push_back(kinValueVecTemp(j));
//                 }

//         }

//         std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
// 	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);

// }


// R CONNECTION: Performance testing function for timing operations to R functions
// Benchmarking utility for evaluating computational performance of matrix operations
arma::fvec testTime(int i, arma::fcolvec & m_bVec){
	arma::fvec vec;
	arma::fvec mvec;
	std::cout << "i is " << i << std::endl;
	clock_t t_0;
	t_0 = clock();
        geno.Get_OneSNP_StdGeno(i, &vec);
	clock_t t_1;
	t_1 = clock();
	std::cout << "t_1-t_0 is " << t_1-t_0 << std::endl;
        float val1 = dot(vec,  m_bVec);
	clock_t t_2;
	t_2 = clock();
	std::cout << "t_2-t_1 is " << t_2-t_1 << std::endl;
        mvec = val1 * (vec);
	clock_t t_3;
	t_3 = clock();
	std::cout << "t_3-t_2 is " << t_3-t_2 << std::endl;
	return(mvec);
}



// R CONNECTION: Sparse matrix operations version 2 to R functions
// General sparse matrix manipulations and transformations for mixed model computations
arma::sp_mat gen_sp_v2(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;

    return result;
}



// R CONNECTION: Sparse linear system solver version 2 to R functions
// Alternative sparse solver implementation for mixed model linear algebra
arma::vec gen_spsolve_v2(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    int r = result.n_rows;
    arma::vec y = arma::linspace<arma::vec>(0, 5, r);	
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    arma::vec x = arma::spsolve( result, y ); 
    	
    return x;
}


// R CONNECTION: R-integrated sparse linear system solver to R functions
// Sparse matrix solver designed for seamless integration with R statistical computing
arma::vec gen_spsolve_inR(const arma::sp_mat& a, arma::vec & y) {
    // sparse x sparse -> sparse
    //arma::sp_mat result1 = result * A;
    arma::vec x = arma::spsolve( a, y );

    return x;
}


// R CONNECTION: Returns diagonal elements of kinship matrix to R functions
// Provides self-relationship values (usually 1) for genomic relationship modeling
arma::fvec get_DiagofKin(){
    //int M = geno.getM();
    int Nnomissing = geno.getNnomissing();
        //cout << "MminMAF=" << MminMAF << endl;
        //cout << "M=" << M << endl; 


    arma::fvec x(Nnomissing);

    if(!(geno.setKinDiagtoOne)){
           x  = (*geno.Get_Diagof_StdGeno());
    	   int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
           x = x/MminMAF; 
    }else{
	   x  = arma::ones<arma::fvec>(Nnomissing);	
    }	
    return(x);
}





//The code below is modified from http://gallery.rcpp.org/articles/parallel-inner-product/
struct stdgenoVectorScalorProduct : public Worker
{
   // source vectors
   arma::fvec & m_bout;
   float  y;
   //unsigned int m_N;
   int jthMarker;

   // constructors
   stdgenoVectorScalorProduct(const int jth, const float y, arma::fvec & prodVec)
      : jthMarker(jth), y(y), m_bout(prodVec) {
        //m_N = geno.getNnomissing();
//      m_bout.zeros(m_N);

  }


   // process just the elements of the range I've been asked to

        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                geno.Get_OneSNP_StdGeno(jthMarker, &vec);
                for(unsigned int i = begin; i < end; i++){
                        m_bout[i] = m_bout[i]+vec[i] * y;
                }
        }



};



void getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec) {


   stdgenoVectorScalorProduct stdgenoVectorScalorProduct(jth, y, prodVec);

   unsigned int m_N = geno.getNnomissing();

   parallelFor(0, m_N, stdgenoVectorScalorProduct);

   // return the computed product
}





struct getP_mailman : public Worker
{
        // source vectors
        unsigned int ithMarker;
	unsigned int powVal;
	arma::ivec ithGeno;
        // destination vector
        arma::ivec Psubvec;


        // constructors
        getP_mailman(unsigned int ith, unsigned int mmchunksize)
                : ithMarker(ith){
		ithGeno = Get_OneSNP_Geno(ith);			
		//unsigned int k  = pow(3, mmchunksize);
		unsigned m_N = geno.getNnomissing();
		Psubvec.zeros(m_N);
		unsigned int powNumber = mmchunksize - 1 - ith % mmchunksize; 		
		powVal = pow(3, powNumber);
        }


	// take the square root of the range of elements requested
     void operator()(std::size_t begin, std::size_t end) {

	for(unsigned int j = begin; j < end; j++){
		Psubvec[j] = ithGeno[j] * powVal;
        }		 
     }

};


int computePindex(arma::ivec &ithGeno){
	int a = ithGeno.n_elem;
	int q = 0;
	int baseNum;
	for(unsigned int i = 0; i < a; i++){
		baseNum = pow(3, a - i - 1);
		q = q + ithGeno[i] * baseNum;
	}
	return(q);
}


struct getP_mailman_NbyM : public Worker
{
        // source vectors
        unsigned int jthChunk;
        unsigned int mmchunksize;
        //arma::ivec ithGeno;
        // destination vector
        arma::ivec Psubvec;


        // constructors
        getP_mailman_NbyM(unsigned int jthChunk,unsigned int mmchunksize)
                : jthChunk(jthChunk), mmchunksize(mmchunksize){
                //ithGeno = Get_OneSNP_Geno(ith);
                //unsigned int k  = pow(3, mmchunksize);
                unsigned m_M = geno.getM();
                Psubvec.zeros(m_M);
                //powNumber = mmchunksize - 1 - ith % mmchunksize;
        }


        // take the square root of the range of elements requested
     void operator()(std::size_t begin, std::size_t end) {
	arma::ivec ithGeno;	
	arma::ivec ithGenosub;
	unsigned int jthIndvStart = jthChunk * mmchunksize;
	unsigned int jthIndvEnd = (jthChunk+1) * mmchunksize - 1;
	arma::uvec indvIndex = arma::linspace<arma::uvec>(jthIndvStart, jthIndvEnd);
        for(unsigned int i = begin; i < end; i++){
		ithGeno = Get_OneSNP_Geno(i);		
		ithGenosub = ithGeno.elem(indvIndex);
		Psubvec[i] = computePindex(ithGenosub);
        }
     }

};



// 
//arma::ivec parallelmmGetP(unsigned int ith, unsigned int mmchunksize) {
  
//  	int M = geno.getM();
//	int N = geno.getNnomissing();	
//  	Pvec.zeros(N);

//  	getP_mailman getP_mailman(ith, mmchunksize);
  
//  	parallelFor(0, N, getP_mailman);
 	
//  	return getP_mailman.Psubvec;
//}



void sumPz(arma::fvec & Pbvec, arma::fvec & Ubvec, unsigned int mmchunksize){

        for (int i = 0; i < Pbvec.n_elem; i++){
                std::cout << "i: " << i << " " << Pbvec[i] << std::endl;
        }

        unsigned int d = Pbvec.n_elem;;
        Ubvec.zeros(mmchunksize);
        unsigned int i = 0;
        arma::fvec z0;
        arma::fvec z1;
        arma::fvec z2;
        z0.zeros(d/3);
        z1.zeros(d/3);
        z2.zeros(d/3);

        while(i < mmchunksize){
                d = d / 3;
//              std::cout << "d: " << d << std::endl;
                z0.resize(d);
                z1.resize(d);
                z2.resize(d);

//              arma::uvec indexvec = arma::linspace<arma::uvec>(0, d-1);
                z0 = Pbvec.subvec(0, d-1);
/*
                 for (int j = 0; j < z0.n_elem; j++){
                std::cout << "j: " << j << " " << z0[j] << std::endl;
        }
*/
                //indexvec = arma::linspace<arma::uvec>(d, 2*d-1);
                //z1 = Pbvec.elem(indexvec);
                z1 = Pbvec.subvec(d, 2*d-1);
                //indexvec = arma::linspace<arma::uvec>(2*d, 3*d-1);
                //z2 = Pbvec.elem(indexvec);
                z2 = Pbvec.subvec(2*d, 3*d-1);

                Pbvec.resize(d);
                Pbvec = z0 + z1 + z2;
                Ubvec[i] = sum(z1) + 2*sum(z2);
                i = i + 1;
              std::cout << "i: " << i << std::endl;
              std::cout << "Ubvec[i]: " << Ubvec[i] << std::endl;

        }
}




// R CONNECTION: Mailman algorithm M-by-N matrix-vector multiplication to R functions
// High-performance genotype matrix operations using chunked memory-efficient computation
void mmGetPb_MbyN(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec, arma::fvec & kinbvec) {
	std::cout << "OKKK" << std::endl;
        int M = geno.getM();
        int N = geno.getNnomissing();
	int k = pow(3,mmchunksize);
        arma::ivec Pvec;
	Pvec.zeros(N);
	Pbvec.zeros(k);
	arma::ivec ithGeno;
	ithGeno.ones(N);
	unsigned int Ptemp;
	Ptemp = 1;
	int indL = cthchunk*mmchunksize;
	int indH = (cthchunk+1)*mmchunksize - 1;
	unsigned int j0 = 0;
	//arma::fmat stdGenoMat(mmchunksize, N);
	float ithfreq; 
	float ithinvstd; 
	arma::fvec chunkfreq = geno.alleleFreqVec.subvec(indL, indH);
	arma::fvec chunkinvstd = geno.invstdvVec.subvec(indL, indH); 
	arma::fvec chunkbvec = bvec.subvec(indL, indH); 

	for (int i = indH; i >= indL; i--){
		ithGeno = Get_OneSNP_Geno(i);
		cout << "Ptemp: " << Ptemp << endl;
		//ithfreq = geno.alleleFreqVec(i);
		//ithinvstd = geno.invstdvVec(i);
		Pvec = Pvec + Ptemp * ithGeno; 
		Ptemp = Ptemp * 3;
		//stdGenoMat.row(j) = ithGeno*ithinvstd - 2*ithfreq*ithinvstd;
		//j0 = j0 + 1;

                //unsigned int k  = pow(3, mmchunksize);
                //unsigned m_N = geno.getNnomissing();
                //Psubvec.zeros(m_N);
                //unsigned int powNumber = mmchunksize - 1 - ith % mmchunksize;

	
	//	getP_mailman getP_mailman(i, mmchunksize);
	//	parallelFor(0, N, getP_mailman);
	//	Pvec = Pvec + getP_mailman.Psubvec;
	//	getP_mailman.Psubvec.clear();
  	}
	

	for (int i = 0; i < N; i++){	
//		std::cout << "i: " << i << " " << Pvec[i] << std::endl;	
		Pbvec[Pvec[i]] = Pbvec[Pvec[i]] + bvec[i];
//		std::cout << "Pbvec[Pvec[i]] " << Pbvec[Pvec[i]] << std::endl;
	}

	arma::fvec Gbvectemp;
	sumPz(Pbvec, Gbvectemp, mmchunksize);
	arma::fvec crossKinVec;
	arma::fvec GbvecInvStd = Gbvectemp % chunkinvstd;
        arma::fvec secondTerm = 2*chunkfreq % chunkinvstd * sum(chunkbvec);
        crossKinVec  = GbvecInvStd - secondTerm;

	//getstdgenoVectorScalorProduct(j, crossKinVec[j], kinbvec);
	j0 = 0;
	arma::fvec stdvec;
	for (int i = indL; i <= indH; i++){
		geno.Get_OneSNP_StdGeno(i, &stdvec);
                kinbvec = kinbvec + crossKinVec[j0]*(stdvec);
		j0 = j0 + 1;
	}

//	for (int i = 0; i < k; i++){
//                std::cout << "Pbvec[i]: " << i << " " << Pbvec[i] << std::endl;
//        }

        //return Pbvec;
}


// R CONNECTION: Mailman algorithm N-by-M matrix-vector multiplication to R functions
// Optimized transposed genotype matrix operations for kinship and association analysis
void mmGetPb_NbyM(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec) {

        int M = geno.getM();
        int N = geno.getNnomissing();
        int k = pow(3,mmchunksize);
        arma::ivec Pvec;
        Pvec.zeros(M);
        Pbvec.zeros(k);
	getP_mailman_NbyM getP_mailman_NbyM(cthchunk,mmchunksize);
	parallelFor(0, M, getP_mailman_NbyM);
	Pvec = getP_mailman_NbyM.Psubvec;
	for (int i = 0; i < M; i++){
		Pbvec[Pvec[i]] = Pbvec[Pvec[i]] + bvec[i];
	}
}




// R CONNECTION: Core Mailman matrix multiplication for genotype data to R functions
// Fast matrix-vector products using memory-efficient Mailman algorithm for large-scale GWAS
void muliplyMailman(arma::fvec & bvec, arma::fvec & Gbvec, arma::fvec & kinbvec){
	int M = geno.getM();
        int N = geno.getNnomissing();

        Gbvec.zeros(M);
	std::cout << "Gbvec.n_elem " << Gbvec.n_elem << std::endl;
        unsigned int mmchunksize = ceil(log(N)/log(3));
	std::cout << "mmchunksize " << mmchunksize << std::endl;

        int numchunk = M / mmchunksize; 
	std::cout << "numchunk " << numchunk << std::endl;
        int reschunk = M % mmchunksize;
	std::cout << "reschunk " << reschunk << std::endl;
	//unsigned int indL;
	//unsigned int indH;
	//mmGetPb(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec)
	arma::fvec Pbvec;
	//arma::fvec Gbvectemp;	


	
	//for (unsigned int j = 0; j < 1; j++){
	for (unsigned int j = 0; j < numchunk; j++){
//		std::cout << "j: " << j << std::endl;
		//Pbvec.zeros(M);
		//indL = j*mmchunksize;
		//indH = (j+1)*mmchunksize-1;
//		if(j == 0){
		double wall0ain = get_wall_time();
 		double cpu0ain  = get_cpu_time();
//		}
//		mmGetPb_MbyN(j, mmchunksize, bvec, Pbvec);

//		if(j == 0){

		mmGetPb_MbyN(j, mmchunksize, bvec, Pbvec, kinbvec);



	double wall1ain = get_wall_time();
 double cpu1ain  = get_cpu_time();
 cout << "Wall Time in mmGetPb_MbyN = " << wall1ain - wall0ain << endl;
 cout << "CPU Time  in mmGetPb_MbyN = " << cpu1ain - cpu0ain  << endl;


//}

//		sumPz(Pbvec, Gbvectemp, mmchunksize);

//if(j == 0){
cout << "ith chunk " << j << endl;
//}
		//getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec)

//		Gbvec.subvec(j*mmchunksize, (j+1)*mmchunksize-1) = Gbvectemp;
  	}

        if(reschunk > 0){
			arma::fvec vec;
		//arma::uvec indexvec = arma::linspace<arma::uvec>(M-reschunk, M-1);
		for (unsigned int j = M-reschunk; j < M; j++){
     		           geno.Get_OneSNP_StdGeno(j, &vec);
			kinbvec = kinbvec + arma::dot(vec, bvec) * vec;
		}	
        }

	kinbvec = kinbvec / M;
}



// R CONNECTION: Mailman N-by-M multiplication for transposed operations to R functions
// Efficient transposed genotype matrix multiplication for kinship matrix construction
void muliplyMailman_NbyM(arma::fvec & bvec, arma::fvec & tGbvec){
        int M = geno.getM();
        int N = geno.getNnomissing();

        tGbvec.zeros(N);

        unsigned int mmchunksize = ceil(log(M)/log(3));

        int numchunk = N / mmchunksize;
        int reschunk = N % mmchunksize;
        unsigned int indL;
        unsigned int indH;
        //mmGetPb(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec)
        arma::fvec Pbvec;
        Pbvec.zeros(M);
	arma::fvec tGbvectemp;

        for (unsigned int j = 0; j < numchunk; j++){
                indL = j*mmchunksize;
                indH = (j+1)*mmchunksize-1;
		mmGetPb_NbyM(j, mmchunksize, bvec, Pbvec);
           	sumPz(Pbvec, tGbvectemp, mmchunksize);
                tGbvec.subvec(j*mmchunksize, (j+1)*mmchunksize-1) = tGbvectemp;     
        }

        if(reschunk > 0){
		arma::imat A(reschunk,M);
		A.zeros();
		arma::ivec Gtemp(N);
		Gtemp.zeros();
		arma::ivec Gtemp2(reschunk);
		Gtemp2.zeros();
		arma::uvec indexvec = arma::linspace<arma::uvec>(M - reschunk -1, M);
                for (unsigned int j = 0; j < M; j++){
			Gtemp = Get_OneSNP_Geno(j);
			Gtemp2 = Gtemp.elem(indexvec);
			A.col(j) = Gtemp2;
                }

		Pbvec.elem(indexvec) = Gtemp2 * (bvec.elem(indexvec));
        }
}


// R CONNECTION: Computes frequency over standard deviation vector to R functions
// Calculates normalized allele frequency statistics for genotype standardization
void freqOverStd(arma::fcolvec& freqOverStdVec){
	freqOverStdVec = 2 * (geno.alleleFreqVec) % (geno.invstdvVec);

	 //int M = geno.getM();
/*
	for (unsigned int j = 0; j < M; j++){
		std::cout << "geno.alleleFreqVec " << j << " " << geno.alleleFreqVec[j] << std::endl; 
		std::cout << "geno.invstdvVec " << j << " " << geno.invstdvVec[j] << std::endl; 
		std::cout << "freqOverStdVec " << j << " " << freqOverStdVec[j] << std::endl; 
               }
*/

}

// R CONNECTION: Mailman-based cross-product with kinship matrix to R functions
// Combines genotype matrix operations with kinship relationships using Mailman algorithm
arma::fvec getCrossprodMatAndKin_mailman(arma::fcolvec& bVec){
	std::cout << "b0: " << std::endl;
	int M = geno.getM();
        int N = geno.getNnomissing();
	arma::fvec Gbvec;


	double wall0in = get_wall_time();
 	double cpu0in  = get_cpu_time();
 	arma::fvec kinbvec;
        kinbvec.zeros(N);

	muliplyMailman(bVec, Gbvec, kinbvec);


double wall1in = get_wall_time();
 double cpu1in  = get_cpu_time();
 cout << "Wall Time in muliplyMailman = " << wall1in - wall0in << endl;
 cout << "CPU Time  in muliplyMailman = " << cpu1in - cpu0in  << endl;



//	for (unsigned int j = 0; j < M; j++){
//                std::cout << "Gbvec " << j << " " << Gbvec[j] << std::endl;
//               }
/*
//	std::cout << "b: " << std::endl;
	arma::fvec freqOverStdVec;
//	std::cout << "a: " << std::endl;
	freqOverStd(freqOverStdVec);
//	std::cout << "c: " << std::endl;
	arma::fvec crossKinVec;
	arma::fvec GbvecInvStd = Gbvec % (geno.invstdvVec);
	arma::fvec secondTerm = freqOverStdVec * sum(bVec);
	crossKinVec  = GbvecInvStd - secondTerm;

double wall2in = get_wall_time();
 double cpu2in  = get_cpu_time();
 cout << "Wall Time in Gtb = " << wall2in - wall1in << endl;
 cout << "CPU Time  in Gtb = " << cpu2in - cpu1in  << endl;


	 for (unsigned int j = 0; j < M; j++){
                std::cout << "GbvecInvStd " << j << " " << GbvecInvStd[j] << std::endl;
                std::cout << "secondTerm " << j << " " << secondTerm[j] << std::endl;
		std::cout << "crossKinVec " << j << " " << crossKinVec[j] << std::endl;
               }
*/
/*	
	arma::fvec kinbvec;
	kinbvec.zeros(N);

	for (unsigned int j = 0; j < M; j++){
		getstdgenoVectorScalorProduct(j, crossKinVec[j], kinbvec);
	}


double wall3in = get_wall_time();
 double cpu3in  = get_cpu_time();
 cout << "Wall Time in getstdgenoVectorScalorProduct = " << wall3in - wall2in << endl;
 cout << "CPU Time  in getstdgenoVectorScalorProduct = " << cpu3in - cpu2in  << endl;



	kinbvec = kinbvec / M;
*/	
        return(kinbvec);

}


// R CONNECTION: Returns diagonal of genomic relationship matrix to R functions
// Provides self-relationship values from GRM for mixed model variance component estimation
arma::fvec get_GRMdiagVec(){
  int mMarker = gettotalMarker(); 
  int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        //cout << "MminMAF=" << MminMAF << endl;

  arma::fvec diagGRMVec = (*geno.Get_Diagof_StdGeno())/MminMAF;
  return(diagGRMVec);
}


// R CONNECTION: Sets minimum allele frequency threshold for GRM construction from R functions
// Configures MAF cutoff for including markers in genomic relationship matrix
void setminMAFforGRM(float minMAFforGRM){
  minMAFtoConstructGRM = minMAFforGRM;
}


// R CONNECTION: Sets maximum missing rate threshold for GRM markers from R functions
// Configures quality control threshold for marker inclusion in kinship analysis
void setmaxMissingRateforGRM(float maxMissingforGRM){
  geno.maxMissingRate = maxMissingforGRM;
}



// R CONNECTION: Configures diagonal standardized genotype matrix for LOCO from R functions
// Sets up leave-one-chromosome-out standardized genotype diagonal elements
void set_Diagof_StdGeno_LOCO(){

      
  int Nnomissing = geno.getNnomissing();
  int chrlength = geno.startIndexVec.n_elem;
  (geno.mtx_DiagStd_LOCO).zeros(Nnomissing, chrlength);
  (geno.Msub_MAFge_minMAFtoConstructGRM_byChr).zeros(chrlength);
//  std::cout << "debug1" << std::endl;
    int starti, endi;
    arma::fvec * temp = &geno.m_OneSNP_StdGeno;
for(size_t k=0; k< chrlength; k++){
   starti = geno.startIndexVec[k];
   endi = geno.endIndexVec[k];
//  std::cout << "debug2" << std::endl;
  if((starti != -1) && (endi != -1)){
  	for(int i=starti; i<= endi; i++){
         		geno.Get_OneSNP_StdGeno(i, temp);
	 		(geno.mtx_DiagStd_LOCO).col(k) = (geno.mtx_DiagStd_LOCO).col(k) + (*temp) % (*temp);
	 		geno.Msub_MAFge_minMAFtoConstructGRM_byChr[k] = geno.Msub_MAFge_minMAFtoConstructGRM_byChr[k] + 1;

  	}
  (geno.mtx_DiagStd_LOCO).col(k) = *geno.Get_Diagof_StdGeno() -  (geno.mtx_DiagStd_LOCO).col(k);
  }
}	
}

/*

// R CONNECTION: Sets MAC thresholds for variance ratio categories from R functions
// Configures minor allele count ranges for variance component ratio estimation
void setminMAC_VarianceRatio(arma::fvec  t_cateVarRatioMinMACVecExclude, arma::fvec  t_cateVarRatioMaxMACVecInclude){
  g_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
  g_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
}
*/



void setminMAC_VarianceRatio(float t_minMACVarRatio, float t_maxMACVarRatio, bool t_isVarianceRatioinGeno){ 
	geno.g_minMACVarRatio = t_minMACVarRatio;
	geno.g_maxMACVarRatio = t_maxMACVarRatio;
	geno.isVarRatio = t_isVarianceRatioinGeno;
	std::cout << "geno.g_minMACVarRatio " << geno.g_minMACVarRatio << std::endl;
	std::cout << "geno.g_maxMACVarRatio " << geno.g_maxMACVarRatio << std::endl;	
}

//  
//int getNumofMarkersforGRM(){
//  int a = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
//  return(a);
//}


// SURVIVAL ANALYSIS FUNCTIONS


// Rcpp::List GetIndexofCases(const arma::vec& status, const arma::vec& time) {
//     int n = time.n_elem;
    
//     // Create data frame equivalent structure
//     arma::vec timeVec = time;
//     arma::uvec orgIndex = arma::linspace<arma::uvec>(0, n-1, n);
//     arma::vec statusVec = status;
    
//     // Sort by time
//     arma::uvec sortIdx = arma::sort_index(timeVec);
//     timeVec = timeVec(sortIdx);
//     statusVec = statusVec(sortIdx);
//     arma::uvec orgIndexSorted = orgIndex(sortIdx);
    
//     // Create newIndex (0-based indexing)
//     arma::uvec newIndex = arma::linspace<arma::uvec>(0, n-1, n);
    
//     // Find case indices (status == 1)
//     arma::uvec caseIndex = arma::find(statusVec == 1);
//     arma::uvec caseIndexwithTies = caseIndex;
    
//     // Handle ties
//     for(arma::uword i = 1; i < caseIndex.n_elem; i++) {
//         if(timeVec(caseIndex(i)) == timeVec(caseIndex(i-1))) {
//             caseIndexwithTies(i) = caseIndexwithTies(i-1);
//         }
//     }
    
//     // Get unique time indices
//     arma::uvec uniqTimeIndex = arma::unique(caseIndexwithTies);
//     arma::vec uniqTimeVec = timeVec(uniqTimeIndex);
    
//     // Create newIndexWithTies
//     arma::uvec newIndexWithTies = newIndex;
//     for(int i = 0; i < n; i++) {
//         for(arma::uword j = 0; j < uniqTimeVec.n_elem; j++) {
//             if(timeVec(i) == uniqTimeVec(j)) {
//                 newIndexWithTies(i) = uniqTimeIndex(j);
//                 break;
//             }
//         }
//     }
    
//     return Rcpp::List::create(
//         Rcpp::Named("timedata") = Rcpp::DataFrame::create(
//             Rcpp::Named("time") = timeVec,
//             Rcpp::Named("orgIndex") = orgIndexSorted,
//             Rcpp::Named("status") = statusVec,
//             Rcpp::Named("newIndex") = newIndex,
//             Rcpp::Named("newIndexWithTies") = newIndexWithTies
//         ),
//         Rcpp::Named("caseIndex") = caseIndex,
//         Rcpp::Named("caseIndexwithTies") = caseIndexwithTies,
//         Rcpp::Named("uniqTimeIndex") = uniqTimeIndex
//     );
// }


arma::vec GetdenominN(const arma::uvec& uniqTimeIndex, 
                      const arma::vec& lin_pred_new, 
                      const arma::uvec& newIndexWithTies, 
                      const arma::uvec& caseIndexwithTies, 
                      const arma::uvec& orgIndex) {
    
    arma::vec explin = arma::exp(lin_pred_new);
    arma::vec demonVec(uniqTimeIndex.n_elem);
    
    for(arma::uword i = 0; i < uniqTimeIndex.n_elem; i++) {
        int nc = explin.n_elem;
        int ntie = arma::sum(caseIndexwithTies == uniqTimeIndex(i));
        
        // Find indices where newIndexWithTies >= uniqTimeIndex[i]
        arma::uvec riskSet = arma::find(newIndexWithTies >= uniqTimeIndex(i));
        
        double denomSum = 0.0;
        for(arma::uword j = 0; j < riskSet.n_elem; j++) {
            denomSum += explin(orgIndex(riskSet(j)));
        }
        
        demonVec(i) = ntie / (denomSum * denomSum);
    }
    
    return demonVec;
}


arma::vec GetdenominLambda0(const arma::uvec& caseIndexwithTies, 
                            const arma::vec& lin_pred_new, 
                            const arma::uvec& newIndexWithTies) {
    
    arma::vec explin = arma::exp(lin_pred_new);
    arma::vec demonVec(caseIndexwithTies.n_elem);
    
    for(arma::uword i = 0; i < caseIndexwithTies.n_elem; i++) {
        arma::uvec riskSet = arma::find(newIndexWithTies >= caseIndexwithTies(i));
        double denomSum = arma::sum(explin(riskSet));
        demonVec(i) = 1.0 / denomSum;
    }
    
    return demonVec;
}


arma::vec GetLambda0(const arma::vec& lin_pred, const Rcpp::List& inC) {
    Rcpp::DataFrame timedata = Rcpp::as<Rcpp::DataFrame>(inC["timedata"]);
    arma::uvec orgIndex = timedata["orgIndex"];
    arma::uvec caseIndexwithTies = inC["caseIndexwithTies"];
    arma::uvec newIndexWithTies = timedata["newIndexWithTies"];
    
    arma::vec lin_pred_new(lin_pred.n_elem);
    for(arma::uword i = 0; i < orgIndex.n_elem; i++) {
        lin_pred_new(i) = lin_pred(orgIndex(i));
    }
    
    arma::vec demonVec = GetdenominLambda0(caseIndexwithTies, lin_pred_new, newIndexWithTies);
    arma::vec Lambda0_vec = arma::cumsum(demonVec);
    
    return Lambda0_vec;
}

// COVARIATE TRANSFORMATION FUNCTIONS


// Rcpp::List Covariate_Transform(arma::mat& X, double tol = 1e-7) {
//     int n = X.n_rows;
//     int p = X.n_cols;
    
//     // Check for multicollinearity using QR decomposition
//     arma::mat Q, R;
//     arma::qr_econ(Q, R, X);
    
//     // Find columns to keep based on diagonal elements of R
//     arma::uvec keep_cols;
//     for(int j = 0; j < p; j++) {
//         if(std::abs(R(j, j)) > tol) {
//             keep_cols.insert_rows(keep_cols.n_elem, 1);
//             keep_cols(keep_cols.n_elem - 1) = j;
//         }
//     }
    
//     // Extract the relevant columns and QR components
//     arma::mat X_reduced = X.cols(keep_cols);
//     arma::mat Q_final, R_final;
//     arma::qr_econ(Q_final, R_final, X_reduced);
    
//     // Transform the design matrix
//     arma::mat X_transformed = Q_final * arma::sqrt(arma::eye(Q_final.n_cols, Q_final.n_cols) * n);
    
//     return Rcpp::List::create(
//         Rcpp::Named("X_transformed") = X_transformed,
//         Rcpp::Named("Q") = Q_final,
//         Rcpp::Named("R") = R_final,
//         Rcpp::Named("keep_cols") = keep_cols,
//         Rcpp::Named("rank") = keep_cols.n_elem
//     );
// }


arma::vec Covariate_Transform_Back(const arma::vec& coeff_transformed, 
                                   const arma::mat& Q, 
                                   const arma::mat& R,
                                   const arma::uvec& keep_cols,
                                   int original_p) {
    
    // Back-transform coefficients
    arma::vec coeff_reduced = arma::solve(arma::trimatu(R), Q.t() * coeff_transformed);
    
    // Expand to original dimension
    arma::vec coeff_original = arma::zeros(original_p);
    for(arma::uword i = 0; i < keep_cols.n_elem; i++) {
        coeff_original(keep_cols(i)) = coeff_reduced(i);
    }
    
    return coeff_original;
}

// PCG SOLVER FUNCTIONS


arma::vec pcg(const arma::mat& A, const arma::vec& b, const arma::vec& M_inv, 
              double tol = 1e-6, int maxiter = 1000) {
    
    int n = A.n_rows;
    arma::vec x = arma::zeros(n);
    arma::vec r = b - A * x;
    arma::vec z = M_inv % r;  // Element-wise multiplication for diagonal preconditioning
    arma::vec p = z;
    
    double rsold = arma::dot(r, z);
    
    for(int iter = 0; iter < maxiter; iter++) {
        arma::vec Ap = A * p;
        double alpha = rsold / arma::dot(p, Ap);
        
        x += alpha * p;
        r -= alpha * Ap;
        
        double rnorm = arma::norm(r, 2);
        if(rnorm < tol) {
            break;
        }
        
        z = M_inv % r;
        double rsnew = arma::dot(r, z);
        double beta = rsnew / rsold;
        
        p = z + beta * p;
        rsold = rsnew;
    }
    
    return x;
}


arma::vec pcgSparse(const arma::sp_mat& A, const arma::vec& b, const arma::vec& M_inv,
                   double tol = 1e-6, int maxiter = 1000) {
    
    int n = A.n_rows;
    arma::vec x = arma::zeros(n);
    arma::vec r = b - A * x;
    arma::vec z = M_inv % r;
    arma::vec p = z;
    
    double rsold = arma::dot(r, z);
    
    for(int iter = 0; iter < maxiter; iter++) {
        arma::vec Ap = A * p;
        double alpha = rsold / arma::dot(p, Ap);
        
        x += alpha * p;
        r -= alpha * Ap;
        
        double rnorm = arma::norm(r, 2);
        if(rnorm < tol) {
            break;
        }
        
        z = M_inv % r;
        double rsnew = arma::dot(r, z);
        double beta = rsnew / rsold;
        
        p = z + beta * p;
        rsold = rsnew;
    }
    
    return x;
}

// COEFFICIENT ESTIMATION FUNCTIONS


arma::vec Get_Coef(arma::mat& X, arma::vec& y, arma::vec& mu, arma::vec& mu2, 
                   arma::vec& coeffs, std::string traitType = "binary",
                   double tol = 1e-6, int maxiter = 30) {
    
    int n = X.n_rows;
    int p = X.n_cols;
    
    for(int iter = 0; iter < maxiter; iter++) {
        arma::vec eta = X * coeffs;
        arma::vec var_mu(n);
        arma::vec dmu_deta(n);
        
        if(traitType == "binary") {
            // Logistic regression
            mu = 1.0 / (1.0 + arma::exp(-eta));
            // Prevent numerical issues
            mu = arma::clamp(mu, 1e-8, 1 - 1e-8);
            var_mu = mu % (1.0 - mu);
            dmu_deta = var_mu;
        } else if(traitType == "quantitative") {
            // Linear regression
            mu = eta;
            var_mu.fill(1.0);
            dmu_deta.fill(1.0);
        }
        
        // Working weights and working response
        arma::vec W = (dmu_deta % dmu_deta) / var_mu;
        arma::vec working_y = eta + (y - mu) / dmu_deta;
        
        // Weighted least squares update
        arma::mat XtWX = X.t() * arma::diagmat(W) * X;
        arma::vec XtWz = X.t() * (W % working_y);
        
        arma::vec delta_coeff = arma::solve(XtWX, XtWz) - coeffs;
        coeffs += delta_coeff;
        
        // Check convergence
        if(arma::norm(delta_coeff) < tol) {
            break;
        }
    }
    
    // Update mu2 for variance calculation if needed
    if(traitType == "binary") {
        arma::vec eta = X * coeffs;
        mu = 1.0 / (1.0 + arma::exp(-eta));
        mu = arma::clamp(mu, 1e-8, 1 - 1e-8);
        mu2 = mu % (1.0 - mu);
    } else {
        mu2.fill(1.0);
    }
    
    return coeffs;
}

// SCORE TEST FUNCTIONS


Rcpp::List ScoreTest_NULL_Model(const arma::mat& X, const arma::vec& y, 
                                const arma::vec& mu, const arma::vec& mu2,
                                const arma::mat& Sigma_i, const arma::mat& Sigma_iX) {
    
    int n = X.n_rows;
    int p = X.n_cols;
    
    // Compute P1 matrix: Sigma_i - Sigma_iX (X^T Sigma_i X)^(-1) X^T Sigma_i
    arma::mat XtSigma_iX = X.t() * Sigma_i * X;
    arma::mat XtSigma_iX_inv = arma::inv_sympd(XtSigma_iX);
    arma::mat P1 = Sigma_i - Sigma_iX * XtSigma_iX_inv * Sigma_iX.t();
    
    // Compute residuals
    arma::vec res = y - mu;
    
    // Compute variance matrix components
    arma::mat P2 = P1;
    if(mu2.n_elem > 0) {
        // For non-identity variance (e.g., binary traits)
        P2 = arma::diagmat(arma::sqrt(mu2)) * P1 * arma::diagmat(arma::sqrt(mu2));
    }
    
    return Rcpp::List::create(
        Rcpp::Named("P1") = P1,
        Rcpp::Named("P2") = P2,
        Rcpp::Named("residuals") = res,
        Rcpp::Named("mu") = mu,
        Rcpp::Named("mu2") = mu2
    );
}


Rcpp::List ScoreTest_NULL_Model_survival(const arma::mat& X, const arma::vec& y,
                                         const arma::vec& time, const arma::vec& status,
                                         const arma::vec& lin_pred, const Rcpp::List& inC) {
    
    // This is a simplified version - full Cox model implementation would be more complex
    int n = X.n_rows;
    
    // Compute score components for survival model
    arma::vec Lambda0 = GetLambda0(lin_pred, inC);
    
    // Compute martingale residuals (simplified)
    arma::vec mart_res = status - Lambda0;
    
    // Information matrix (simplified)
    arma::mat Info = X.t() * X;  // Simplified - should be Fisher information
    
    return Rcpp::List::create(
        Rcpp::Named("martingale_residuals") = mart_res,
        Rcpp::Named("information_matrix") = Info,
        Rcpp::Named("Lambda0") = Lambda0
    );
}
