//ExtractAcceptance.h 
// By Jixie Zhang, 

#ifndef _ExtractAcceptance_
#define _ExtractAcceptance_

#include <string>
#include "ReadSingleArm.h"

#define ExtractAcceptance_Debug 0
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
class ExtractAcceptance:public ReadSingleArm
{
public:
  ExtractAcceptance(const char* filename, int det, int type=1);
  virtual ~ExtractAcceptance();

  void  BeginOfRun();
  void  Run();
  void  EndOfRun();
  ///////////////////////////////////////////////////////////////////////////
  
public:
 
  int  iYtarIdx,iDeltaIdx,iThetaIdx,iPhiIdx;

  int ****N_Inc;   //[kYtarBinNum][kDeltaBinNum][kThetaBinNum][kPhiBinNum];
  int ****N_Inc_true;//[kYtarBinNum][kDeltaBinNum][kThetaBinNum][kPhiBinNum];
  double ****Acc_Inc;//[kYtarBinNum][kDeltaBinNum][kThetaBinNum][kPhiBinNum];
 
private:
  void FillAccTree();
  void Reset();  //reset the tree variables
  void DoAccCal();
  void PrintACCTable(int level, const char *filename, char *headblock, int ****N_det, int ****N_true);
  void PrintACCTable(int level, const char *filename, char *headblock, double ****N_det, double ****N_true);

private:
  int mDet,mType;
  std::string mDetName;
};
#endif //_ExtractAcceptance_

