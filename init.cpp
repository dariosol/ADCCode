
void init(){

  //	gROOT->ProcessLine(".L Risultati.h+");
	gROOT->ProcessLine(".L Risultati.cpp+");
	//gROOT->ProcessLine(".L Header.h+");
	gROOT->ProcessLine(".L Header.cpp+");
	gROOT->ProcessLine(".L istogramma.cpp");
	gROOT->ProcessLine(".L fit.cpp");
	gROOT->ProcessLine(".L fit_original.cpp");
	gROOT->ProcessLine(".L fft.cpp");
	gROOT->ProcessLine(".L fftaverage.cpp");
	gROOT->ProcessLine(".L fft_original.cpp");
	//gROOT->ProcessLine(".L LetturaFile.cpp");
	//	gROOT->ProcessLine(".L GetHeader.cpp");
	gROOT->ProcessLine(".L txtToROOT.cpp");
	//	gROOT->ProcessLine(".L txtToHistROOT.cpp");
	//gROOT->ProcessLine(".L CoprimeFinder.cpp");
	//	gROOT->ProcessLine(".L scriviTxt.cpp");
	gROOT->ProcessLine(".L readeth.C");
	//gROOT->ProcessLine(".L readcsv.C");
	//	gROOT->ProcessLine(".L readeth.C");
	gROOT->ProcessLine(".L calibrazione.cpp");
	gROOT->ProcessLine(".L fftcopy.C");
	gROOT->ProcessLine(".L doallstep.C");
}
