
void r_Signal()
{
//gROOT->ProcessLine(".L RoccoR.cc++");
gSystem->Load("libFWCoreFWLite.so");
//gROOT->ProcessLine(".L BTagCalibrationStandalone.cpp+");
gROOT->ProcessLine(".L Signal.C");
gROOT->ProcessLine("Signal t");
gROOT->ProcessLine("t.Loop()");
gROOT->ProcessLine(".q");
}


