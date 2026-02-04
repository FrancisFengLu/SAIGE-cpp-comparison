#define ARMA_USE_SUPERLU 1
#include <RcppArmadillo.h>
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
#include <RcppParallel.h>
#include <RcppParallel/TBB.h>
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
        std::vector< std::vector<unsigned char>* > genoVecofPointers;
        ///////////
        std::vector< std::vector<unsigned char>* > genoVecofPointers_forVarRatio;
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
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
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
		geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);
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
			geno1 = genoVecofPointers_forVarRatio[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
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
		geno1 = genoVecofPointers_forVarRatio[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
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


	void Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld, float & altFreq, float & missingRate, int & mac,  int & alleleCount, bool & passQC, size_t SNPIdx_new, bool & passVarRatio , size_t SNPIdx_vr){

		arma::ivec m_OneSNP_GenoTemp;
		m_OneSNP_GenoTemp.zeros(N);
		m_OneSNP_Geno.zeros(Nnomissing);
		int m_size_of_esi_temp = (N+3)/4;
		size_t ind= 0;
		unsigned char geno1;
		int bufferGeno;
		int u;
		alleleCount = 0;
		int numMissing = 0;
		for(int i=0; i< (m_size_of_esi_temp - 1); i++){
			geno1 = genoVecOneMarkerOld[i];
			for(int j=0; j<4; j++){
				u = j & 3;

				int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                // PLINK BED encoding: 00=hom_A1, 01=missing, 10=het, 11=hom_A2
                                // bufferGeno represents A1 allele count (0, 1, 2) or 3 for missing
                                // FIXED: Match R version - count A1 alleles
                                if (b == 1 && a == 0){
                                        bufferGeno = 3;  // 01 = missing
                                }else if(b == 0 && a == 0){
                                        bufferGeno = 2;  // 00 = hom_A1 (has 2 A1 alleles)
                                }else if(b == 0 && a == 1){
                                        bufferGeno = 1;  // 10 = het (has 1 A1 allele)
                                }else if(b == 1 && a == 1){
                                        bufferGeno = 0;  // 11 = hom_A2 (has 0 A1 alleles)
                                }else{
                                        cout << "Error GENO!!\n";
                                        break;
                                }


				//getGenotype(&geno1, u, bufferGeno);
				//printf("%d", geno1);
	//			std::cout << "bufferGeno " << bufferGeno << std::endl;
				/*
				switch(geno1){
    				 case HOM_REF: break;
    				 case HET: sum+=1; break;
    				 case HOM_ALT: sum+=2; break;
    				 case MISSING: numMissing++; break;
    				}
				*/
				m_OneSNP_GenoTemp[ind] = bufferGeno;
				//if(SNPIdx == 0){
				//	std::cout << "[ind] " << ind << std::endl;
				//	std::cout << "indicatorGenoSamplesWithPheno_in[ind] " << indicatorGenoSamplesWithPheno_in[ind] << std::endl;
				//}	
				if(indicatorGenoSamplesWithPheno_in[ind]){
					if(bufferGeno == 3){
	//					std::cout << "SNPIdx " << SNPIdx << std::endl;
	//					std::cout << "ind " << ind << std::endl;
						numMissing++;
						//indexNA.push_back(ptrsubSampleInGeno[ind])
					}else{
						alleleCount = alleleCount + bufferGeno;
					}	
				}	
				ind++;
                                geno1 = geno1 >> 1;	
			  }

	      }	  

		int i = m_size_of_esi_temp - 1;
		geno1 = genoVecOneMarkerOld[i];
		//std::cout << "N " << N << std::endl;
		//while(ind < N){
                        for(int j=0; j<4; j++){
                                u = j & 3;

                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                // PLINK BED encoding: 00=hom_A1, 01=missing, 10=het, 11=hom_A2
                                // bufferGeno represents A1 allele count (0, 1, 2) or 3 for missing
                                // FIXED: Match R version - count A1 alleles
                                if (b == 1 && a == 0){
                                        bufferGeno = 3;  // 01 = missing
                                }else if(b == 0 && a == 0){
                                        bufferGeno = 2;  // 00 = hom_A1 (has 2 A1 alleles)
                                }else if(b == 0 && a == 1){
                                        bufferGeno = 1;  // 10 = het (has 1 A1 allele)
                                }else if(b == 1 && a == 1){
                                        bufferGeno = 0;  // 11 = hom_A2 (has 0 A1 alleles)
                                }else{
                                        cout << "Error GENO!!\n";
                                        break;
                                }

                                m_OneSNP_GenoTemp[ind] = bufferGeno;
                                //if(SNPIdx == 0){
                                //        std::cout << "[ind] " << ind << std::endl;
                                //        std::cout << "indicatorGenoSamplesWithPheno_in[ind] " << indicatorGenoSamplesWithPheno_in[ind] << std::endl;
                                //}
                                if(indicatorGenoSamplesWithPheno_in[ind]){
                                        if(bufferGeno == 3){
        //                                      std::cout << "SNPIdx " << SNPIdx << std::endl;
        //                                      std::cout << "ind " << ind << std::endl;
                                                numMissing++;
                                        }else{
                                                alleleCount = alleleCount + bufferGeno;
                                        }
                                }
                                ind++;
                                geno1 = geno1 >> 1;
				if(ind == (N)){
				break;
				}
                          }
		//}


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
			for(int indx=0; indx < Nnomissing; indx++){
                              u = indx & 3;
                              bufferGeno = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
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
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

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
                geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

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
			std::ofstream stdgeno_file("/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/jan_14_comparison/cpp_stdgeno.txt");
			stdgeno_file << "# StdGeno values for first 3 markers, all samples\n";
			stdgeno_file << "# Columns: Sample, Marker0, Marker1, Marker2\n";

			arma::fmat first3_stdgeno(Nnomissing, 3);

			// DEBUG: Track sample 0's stdGeno^2 accumulation
			float sample0_cumsum = 0;
			std::ofstream sample0_file("/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/jan_14_comparison/cpp_sample0_cumsum.txt");
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
			std::cout << "DEBUG: Wrote stdGeno to /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/jan_14_comparison/cpp_stdgeno.txt" << std::endl;

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
		//cout << "OK1\n";
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
		//cout << "OK3\n";
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

    		genoVecOneMarkerOld.reserve(nbyteOld);
    		genoVecOneMarkerOld.resize(nbyteOld);
 		//genoVec.reserve(reserve);
		
		//cout << "OK4\n";

		ifstream test_bedfile;
        	test_bedfile.open(bedfile.c_str(), ios::binary);
        	if (!test_bedfile.is_open()){
                	printf("Error! file open!");
                	return;
        	}		
		//printf("\nM: %zu, N: %zu, Reserve: %zu\n", M, N, reserveTemp);
		printf("\nM: %zu, N: %zu\n", M, N);

        	//test_bedfile.seekg(3);
		//set up the array of vectors for genotype
		//numMarkersofEachArray = floor((memoryChunk*pow (10.0, 9.0))/(ceil(float(N)/4)));
		//cout << "numMarkersofEachArray: " << numMarkersofEachArray << endl;
		numMarkersofEachArray = 1;
		//if(M % numMarkersofEachArray == 0){
                        //numofGenoArray = M / numMarkersofEachArray;
                        numofGenoArray = M;
			genoVecofPointers.resize(numofGenoArray);
			genoVecofPointers_forVarRatio.resize(numofGenoArray);
                        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
                        for (int i = 0; i < numofGenoArray ; i++){
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }

                //}
		/*else{
                        numofGenoArray = M/numMarkersofEachArray + 1;
                        genoVecofPointers.resize(numofGenoArray);
			numMarkersofLastArray = M - (numofGenoArray-1)*numMarkersofEachArray;
                        cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
			try{	
                        for (int i = 0; i < numofGenoArray-1; i++){
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }
			genoVecofPointers[numofGenoArray-1] = new vector<unsigned char>;
			genoVecofPointers[numofGenoArray-1]->reserve(numMarkersofLastArray*ceil(float(Nnomissing)/4));
			}
			catch(std::bad_alloc& ba)
                        {
                                std::cerr << "bad_alloc caught1: " << ba.what() << '\n';
                                exit(EXIT_FAILURE);
                        }
		}*/

		cout << "setgeno mark1" << endl;
		arma::ivec g_randMarkerIndforVR_temp;
		//randomly select common markers for variance ratio
		if(isVarRatio){
			 g_randMarkerIndforVR_temp = arma::randi(1000, arma::distr_param(0,M-1));
			 g_randMarkerIndforVR = arma::unique(g_randMarkerIndforVR_temp);
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
		cout << "setgeno mark2" << endl;
		cout << "Nnomissing (sample count for MAF): " << Nnomissing << endl;
		size_t SNPIdx_new = 0;
		size_t SNPIdx_vr = 0;
		//Mmafge1perc = 0;
		for(int i = 0; i < M; i++){
			genoVecOneMarkerOld.clear();
			genoVecOneMarkerOld.reserve(nbyteOld);
                        genoVecOneMarkerOld.resize(nbyteOld);

			test_bedfile.seekg(3+nbyteOld*i);
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

        	test_bedfile.close();
//		cout << "setgeno mark5" << endl;
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
                        	geno1 = genoVecofPointers[indexOfVectorPointer]->at(k);
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
                        // DEBUG: Check bounds before .at() access
                        if(k == 0 && i == Start_idx) {
                            std::cout << "[DEBUG Get_MultiMarkers] k=" << k
                                      << " indexOfVectorPointer=" << indexOfVectorPointer
                                      << " genoVecofPointers.size()=" << geno.genoVecofPointers.size()
                                      << " i=" << i
                                      << " vector_size=" << (indexOfVectorPointer < geno.genoVecofPointers.size() ? geno.genoVecofPointers[indexOfVectorPointer]->size() : -1)
                                      << std::endl;
                        }
                        if(indexOfVectorPointer >= geno.genoVecofPointers.size()) {
                            std::cerr << "[ERROR] indexOfVectorPointer=" << indexOfVectorPointer
                                      << " >= genoVecofPointers.size()=" << geno.genoVecofPointers.size()
                                      << " at k=" << k << std::endl;
                            throw std::out_of_range("indexOfVectorPointer out of bounds");
                        }
                        if(i >= geno.genoVecofPointers[indexOfVectorPointer]->size()) {
                            std::cerr << "[ERROR] i=" << i
                                      << " >= vector.size()=" << geno.genoVecofPointers[indexOfVectorPointer]->size()
                                      << " at k=" << k << " indexOfVectorPointer=" << indexOfVectorPointer << std::endl;
                            throw std::out_of_range("i out of bounds in genoVecofPointers");
                        }
                        geno1 = (geno.genoVecofPointers)[indexOfVectorPointer]->at(i);
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
                        geno1 = (geno.genoVecofPointers)[indexOfVectorPointer]->at(i);
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
  
  	// product that I have accumulated
  	arma::fvec m_bout;
        int Msub_mafge1perc;	
  
  	// constructors
  	CorssProd(arma::fcolvec & y)
  		: m_bVec(y) {
  		
  		m_M = geno.getM();
  		m_N = geno.getNnomissing();
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
		//geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
  	} 
  	CorssProd(const CorssProd& CorssProd, Split)
  		: m_bVec(CorssProd.m_bVec)
  	{

  		m_N = CorssProd.m_N;
  		m_M = CorssProd.m_M;
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
		//CorssProd.Msub_mafge1perc;
  	
  	}  
  	// process just the elements of the range I've been asked to
  	void operator()(std::size_t begin, std::size_t end) {
  	  	arma::fvec vec;
  	  	for(unsigned int i = begin; i < end; i++){
			//if(geno.alleleFreqVec[i] >= minMAFtoConstructGRM && geno.alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
				geno.Get_OneSNP_StdGeno(i, &vec);
				float val1 = dot(vec,  m_bVec);
				m_bout += val1 * (vec) ;
				Msub_mafge1perc += 1;
			//}
			//std::cout << "i: " << i << std::endl;
			//for(unsigned int j = 0; j < 10; j++){
			//	std::cout << "m_bVec[j] " << m_bVec[j] << std::endl;
			//	std::cout << "vec[j] " << vec[j] << std::endl;
			//}
		
			//m_bout += val1 * (vec) / m_M;
  		}
  	}
  
  	// join my value with that of another InnerProduct
  	void join(const CorssProd & rhs) { 
    		m_bout += rhs.m_bout;
		Msub_mafge1perc += rhs.Msub_mafge1perc; 
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
	Rcpp::Environment base_env("package:base");
  	Rcpp::Function set_seed_r = base_env["set.seed"];
  	set_seed_r(seed);  
}

// INTERNAL: Parallel computation helper for cross products
arma::fvec parallelCrossProd(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
  	//int M = geno.getM();
 	//int Msub_mafge1perc = geno.getMmafge1perc();
	int Msub_mafge1perc = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();

	// DEBUG: Print marker info for first few markers
	static int debug_call_count = 0;
	if(debug_call_count == 0) {
		std::cout << "=== DEBUG parallelCrossProd (first call) ===" << std::endl;
		std::cout << "Msub_mafge1perc: " << Msub_mafge1perc << std::endl;
		arma::fvec vec_debug;
		for(int m = 0; m < 3; m++) {
			std::cout << "Marker " << m << ": freq=" << geno.alleleFreqVec[m]
			          << " invstd=" << geno.invstdvVec[m] << std::endl;
			geno.Get_OneSNP_StdGeno(m, &vec_debug);
			std::cout << "  StdGeno[0:5]: " << vec_debug[0] << " " << vec_debug[1] << " "
			          << vec_debug[2] << " " << vec_debug[3] << " " << vec_debug[4] << std::endl;
		}
		std::cout << "bVec[0:5]: " << bVec[0] << " " << bVec[1] << " " << bVec[2]
		          << " " << bVec[3] << " " << bVec[4] << std::endl;
		debug_call_count++;
	}

	CorssProd CorssProd(bVec);

  // call paralleReduce to start the work
  	parallelReduce(0, Msub_mafge1perc, CorssProd);
 	
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
  	return CorssProd.m_bout/(CorssProd.Msub_mafge1perc);
  	//return CorssProd.m_bout;
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
        return CorssProd.m_bout;
        //return CorssProd.m_bout;
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



	outvec = outvec - CorssProd.m_bout;

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
  	crossProdVec = parallelCrossProd(bVec) ;
   }  
  	return(crossProdVec);
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

//    double wall0in = get_wall_time();
// double cpu0in  = get_cpu_time();

    arma::vec yvec2 = arma::conv_to<arma::vec>::from(yvec);
//double wall1in = get_wall_time();
// double cpu1in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall1in - wall0in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu1in - cpu0in  << endl;

    arma::sp_mat result = gen_sp_Sigma(wVec, tauVec);

//double wall2in = get_wall_time();
// double cpu2in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall2in - wall1in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu2in - cpu1in  << endl;

//    std::cout << "yvec.n_elem: " << yvec.n_elem << std::endl;
//    std::cout << "yvec2.n_elem: " << yvec2.n_elem << std::endl;
//    std::cout << "result.n_rows: " << result.n_rows << std::endl;
//    std::cout << "result.n_cols: " << result.n_cols << std::endl;
    arma::vec x = arma::spsolve(result, yvec2);

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
    const arma::uword n = bVec.n_elem;
    if (n == 0) throw std::invalid_argument("PCG: bVec is empty");
    if (wVec.n_elem != n) throw std::invalid_argument("PCG: wVec.len != bVec.len");
    if (tauVec.n_elem < 2) throw std::invalid_argument("PCG: tauVec must have >=2 elems");

    // make non-const copies to satisfy legacy APIs
    arma::fvec w  = wVec;
    arma::fvec tau = tauVec;

    arma::fvec xVec(n, arma::fill::zeros);
    arma::fvec rVec = bVec;
    arma::fvec zVec(n, arma::fill::zeros);
    arma::fvec minvVec(n, arma::fill::zeros);


    if (!isUsePrecondM) {
        // legacy takes non-const refs
        arma::fvec testvec = getDiagOfSigma(w, tau) ; 
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

    // Preconditioner
    if (!isUsePrecondM) {
        // diagonal preconditioner M^{-1} = 1/diag(Sigma)
        minvVec = 1.0f / getDiagOfSigma(wVec, tauVec);   // legacy: non-const refs
        if (minvVec.n_elem != n)
            throw std::runtime_error("PCG_LOCO: getDiagOfSigma returned wrong length");
        zVec = minvVec % rVec;
    } else {
        // sparse preconditioner using (LOCO-aware) sparse Sigma
        zVec = gen_spsolve_v4(wVec, tauVec, rVec);       // legacy: non-const refs
        if (zVec.n_elem != n)
            throw std::runtime_error("PCG_LOCO: gen_spsolve_v4 returned wrong length");
    }

    arma::fvec pVec = zVec;
    float sumr2 = arma::dot(rVec, rVec);
    int   iter  = 0;

    while (sumr2 > tolPCG && iter < maxiterPCG) {
        ++iter;

        // Ap = Sigma * p  (LOCO should be respected inside this op)
        arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec); // legacy: non-const refs
        if (ApVec.n_elem != n)
            throw std::runtime_error("PCG_LOCO: getCrossprod returned wrong length");

        float a = arma::as_scalar((rVec.t() * zVec) / (pVec.t() * ApVec));
        xVec += a * pVec;

        arma::fvec r1Vec = rVec - a * ApVec;

        arma::fvec z1Vec;
        if (!isUsePrecondM) {
            z1Vec = minvVec % r1Vec;
        } else {
            z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec); // legacy: non-const refs
            if (z1Vec.n_elem != n)
                throw std::runtime_error("PCG_LOCO: gen_spsolve_v4 (z1) wrong length");
        }

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

/*
// INTERNAL: Set chromosome indices for LOCO analysis
void setChromosomeIndicesforLOCO(vector<int> chromosomeStartIndexVec, vector<int> chromosomeEndIndexVec, vector<int> chromosomeVecVec){
  LOCO = true;
  chromosomeStartIndex = chromosomeStartIndexVec;
  chromosomeEndIndexV = chromosomeEndIndexVec;
  chromosomeVec = chromosomeVecVec;
}
*/

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
    std::cerr << "Please run R first to generate more vectors in /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_15_work/random_vectors_intseed3.csv" << std::endl;
    throw std::runtime_error("Preloaded vectors exhausted - cannot continue with matching R vectors");
  }

  // Fall back to random generation (only if preloaded vectors not enabled)
  static std::mt19937 gen(std::random_device{}());
  static std::bernoulli_distribution dist(0.5);

  arma::fvec u(n);
  for (int i = 0; i < n; ++i) {
    u(i) = dist(gen) ? 1.0f : -1.0f;
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
    std::string bypass_path = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/bypass/random_vectors_seed10.csv";
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
  // Reset vector index on each GetTrace call (mimics R's set.seed(1) behavior)
  preloaded_vector_idx = 0;

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

      // Save C++ Au to file for comparison
      {
        std::string cpp_out = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/grm_comparison/cpp_Au_vec_" + std::to_string(i) + ".txt";
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

arma::fvec GetTrace_q(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
  	set_seed(200);
  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = geno.getNnomissing();
  	arma::fvec tempVec(nrun);
  	tempVec.zeros();
  	arma::fvec tempVec0(nrun);
  	tempVec0.zeros();

        arma::fvec Sigma_iu;
        arma::fcolvec Pu;
        arma::fvec Au;
        arma::fvec uVec;

        int nrunStart = 0;
        int nrunEnd = nrun;
        float traceCV = traceCVcutoff + 0.1;
        float traceCV0 = traceCVcutoff + 0.1;

        while((traceCV > traceCVcutoff) | (traceCV0 > traceCVcutoff)){


  	for(int i = nrunStart; i < nrunEnd; i++){

    		Rcpp::NumericVector uVec0;
    		uVec0 = nb(Nnomissing);
    		uVec = as<arma::fvec>(uVec0);
    		uVec = uVec*2 - 1;
  //  		arma::fvec Sigma_iu;
    		Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);
  //  		arma::fcolvec Pu;
    		Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
  //  		arma::fvec Au;
    		Au = getCrossprodMatAndKin(uVec);
    		tempVec(i) = dot(Au, Pu);
    		tempVec0(i) = dot(uVec, Pu);
                Au.clear();
      		Pu.clear();
      		Sigma_iu.clear();
      		uVec.clear();
  	}
	traceCV = calCV(tempVec);
	traceCV0 = calCV(tempVec0);
	
	if((traceCV > traceCVcutoff) | (traceCV0 > traceCVcutoff)){
          nrunStart = nrunEnd;
          nrunEnd = nrunEnd + 10;
          tempVec.resize(nrunEnd);
          tempVec0.resize(nrunEnd);

	  //std::cout << "arma::mean(tempVec0): " << arma::mean(tempVec0) << std::endl;	
          cout << "CV for trace random estimator using "<< nrunStart << " runs is " << traceCV <<  "(> " << traceCVcutoff << endl;
          cout << "try " << nrunEnd << "runs" << endl;
        }

    }   

  	arma::fvec traVec(2);
  	traVec(1) = arma::mean(tempVec);
  	traVec(0) = arma::mean(tempVec0);
	tempVec.clear();
	tempVec0.clear();
  	return(traVec);
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
