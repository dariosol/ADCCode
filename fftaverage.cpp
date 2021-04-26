
#include "Risultati.h"
#include "Header.h"

#if !defined (__CLING__) || defined (__ROOTCLING__)
//		C++
#include <iostream>
#include <string>
#include <cmath>
#include <fftw3.h>
//		ROOT
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TVirtualFFT.h>
#include <TGraph.h>
#endif



int fftaverage(

	std::string nomeFile = "datiFFiT.root",
	bool salva = 0,
	double MSPS = 160,
	bool useCalibratedSamples = 1,
	int split = 0		//0: tutti, 1: dispari, 2: pari

	){

  const int AverageValue=16;
  const int CoherentNSamples=16384;

  //definizioni iniziali parametri ADC
	const double nBit = 12;
	const int nCh = (int) TMath::Power(2, nBit);
	const int chMax = nCh - 1;
	const double vmax = 0.6;
	const double vmin = -0.6;
	const double FSR = vmax - vmin;
	const double LSB = FSR / nCh;

	if (split == 1 || split == 2) MSPS = MSPS / 2;
	if (split==3) MSPS = MSPS / 4.;
	const double fSamp = MSPS * 1000000;
	const double tSamp = 1 / fSamp;
	const bool minimalPlotting = 0;
	const bool grafica = 1;
	const bool output = 1;

	//definizioni varie

	const std::string nomeTGraphDefault = "grSamp";
	std::string nomeTGraph = nomeTGraphDefault;
	if ( split == 1 && useCalibratedSamples == 0 ) nomeTGraph += "Odd";
	if ( split == 2 && useCalibratedSamples == 0 ) nomeTGraph += "Even";
	if ( split == 3 && useCalibratedSamples == 0 ) nomeTGraph += "FourByFour";

	if ( useCalibratedSamples ) nomeTGraph += "Calib";
	std::string sampleType = "";
	if ( !useCalibratedSamples && !split ) sampleType += "All";
	if ( split == 1 ) sampleType += "Odd";
	if ( split == 2 ) sampleType += "Even";
	if ( split == 3 ) sampleType += "FourByFour";

	if ( useCalibratedSamples ) sampleType += "Calibrated";

	const std::string freqKey = "freq[Hz]";
	const int maxHarm = 10;		//numero di armoniche per calcolo THD

	//	<LETTURA FILE>

	TFile tf1( nomeFile.c_str(), "update");

	//controllo esistenza del file
	if ( !tf1.IsOpen() ){

		std::cout << "File \"" << nomeFile.c_str() << "\" not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		return 1;

	}

	//apertura del TGraph
	TGraph * gr1 = (TGraph *) tf1.Get( nomeTGraph.c_str() );

	if ( gr1 == nullptr && useCalibratedSamples ) {

		std::cout << "TGraph \"" << nomeTGraph.c_str() << "\" not found " << std::endl;
		std::cout << "Reverting to \"" << nomeTGraphDefault.c_str() << "\" " << std::endl;
		nomeTGraph = nomeTGraphDefault;
		gr1 = (TGraph *) tf1.Get( nomeTGraph.c_str() );

	}

	if ( gr1 == nullptr ){

		std::cout << "TGraph \"" << nomeTGraph.c_str() << "\" not found " << std::endl;
		std::cout << "Aborting execution" << std::endl;
		tf1.Close();
		return 2;
	}

	//		<LETTURA HEADER>
	Header * headerObj = (Header *) tf1.Get( "Header" );

	if ( headerObj == nullptr ){

		std::cout << "Header not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete gr1;
		tf1.Close();
		return 3;

	}

	std::unordered_map < std::string, double > headerMap = headerObj->GetHeader();
	delete headerObj;

	double freq = 0;

	try {

		freq = headerMap.at( freqKey );

	} catch ( const std::out_of_range& invArg ) {

		std::cout << "Key not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete gr1;
		tf1.Close();
		return 4;

	}

	//		</LETTURA HEADER>

	gr1->SetMarkerStyle(7);

	//long long int nPts = gr1->GetN();
	long long int nPts = CoherentNSamples;

	//altre definizioni
	const double omega = 2 * TMath::Pi() * freq;
	const double T = 2 * TMath::Pi() / omega;

	//	</LETTURA FILE>

	//rette
	TF1	* line1 = new TF1( "line1", "0", 0, nPts - 1 );
	line1->SetLineColor(kOrange - 3);
	line1->SetLineWidth(1);
	line1->SetLineStyle(1);

	TF1	* line2 = new TF1( "line2", std::to_string(nCh).c_str(), 0, nPts - 1 );
	line2->SetLineColor(kOrange - 3);
	line2->SetLineWidth(1);
	line2->SetLineStyle(1);

	//looking for maximum:
	vector<int> points;
	vector<int> samples;
	
	for ( int k = 0; k < gr1->GetN(); k++) {
	  
	  double sampling = 0;
	  double code = 0;
	  
	  gr1->GetPoint( k, sampling, code );
	  points.push_back(code);
	  samples.push_back(sampling);
	  
	}
	
	double max = *max_element(points.begin(), points.end());
	max=max;
	cout<<"Max value: "<<max<<endl;

	TH1D * histCamp[AverageValue];
	for(int i=0;i<AverageValue;++i) {	
	  histCamp[i] = new TH1D(Form("histCamp%d",i), "Sampled signal", nPts, 0 - 0.5, nPts - 1 + 0.5);
	}
	cout<<"Here-1!"<<endl;
	//riempimento manuale istogramma perché il metodo GetHistogram() di TGraph non va
	int n=0;
	for(int i=0;i<AverageValue;++i) {	  
	  for ( int k = n; k < nPts+n; k++) {   
	    double sampling = 0;
	    double code = 0; 
	    gr1->GetPoint( k, sampling, code );
	    double HannWin = 0.5 * (1 - cos(2*3.14159265358979323846*(k-n)/(CoherentNSamples-1)));
	    //	    code= HannWin * code;
	    (histCamp[i]->SetBinContent(sampling + 1 -n, code ));

	  }
	  n += CoherentNSamples;
	}

	//	<FFT>
	
	TH1D * histMag[AverageValue];


	for(int i=0;i<AverageValue;++i) {
	  histMag[i] = 0;
	}
	TVirtualFFT::SetTransform(0);	//controllare cosa faccia
	for(int i=0;i<AverageValue;++i) {
	  histMag[i] = (TH1D *) histCamp[i]->FFT( histMag[i], Form("MAG%d",i) );
	  histMag[i]->GetXaxis()->SetRange( 2, nPts / 2 + 1 );
	  histMag[i]->GetXaxis()->SetRange( 0, nPts + 1 );
	}

	//	</FFT>

	//	<ANALISI>
	TH1D * histMagTot = new TH1D(*histMag[0]);
	for(int i =1;i<AverageValue;++i) {
	  histMagTot->Add(histMag[i]);
	  }


	histMagTot->Scale(1./AverageValue);
	histMagTot->GetXaxis()->SetRange( 20, nPts / 2 + 1 );
	int maxBin = histMagTot->GetMaximumBin();
	histMagTot->GetXaxis()->SetRange( 0, nPts + 1 );
       
	cout<<"max Bin: "<<maxBin<<endl;
	TH1D * histMag2 = new TH1D ( *histMagTot );
	TH1D * histMag3 = new TH1D ( *histMagTot );
	//		<NAD>

	double sqSum = 0;

	//sommatoria per NAD
	//Sullo standard c'è scritto di considerare anche le ampiezze oltre la frequenza di Nyquist
	for ( int bin = 1; bin < histMagTot->GetNbinsX() + 1; bin++ ) {		//da 1 ad nPts - 1, ovvero tra il bin 2 ed il bin nPts

		if ( bin != 1 && bin != maxBin && bin != histMagTot->GetNbinsX() - maxBin + 2 ){	//devi anche togliere n - maxBin

			sqSum += histMagTot->GetBinContent( bin ) * histMagTot->GetBinContent( bin );

		}

	}
	cout<<"sqSum "<<sqSum<<endl;
	cout<<"npts "<< nPts<<endl;
	
	double NAD = TMath::Sqrt( sqSum ) / TMath::Sqrt(  nPts * ( nPts - 3 )  );	//ok
	
	//		</NAD>

	//		<A_RMS>

	double Arms = TMath::Sqrt( 2 * histMagTot->GetBinContent( maxBin ) * histMagTot->GetBinContent( maxBin ) ) / nPts;	//ok**

	//		</A_RMS>

	//		<THD>

	double binArm = 0;
	double sqSumHarm = 0;
	int nIn = nPts * freq / fSamp;	//bin di frequenza del segnale in ingresso

	//somma sulle armoniche del segnale(da 2 a maxHarm, inclusa)
	//ma bisogna contare anche quelle sopra la frequenza di Nyquist?? Sì!
	for (int n = 2; n <= maxHarm; n++){

		binArm = ( n * nIn ) % nPts + 1; //il +1 c'è per come son fatti i bin dei TH1
		sqSumHarm += 2 * histMagTot->GetBinContent( binArm ) * histMagTot->GetBinContent( binArm );	//il 2 c'è per le frequenze negative

	}

	sqSumHarm = sqSumHarm / (nPts * nPts);

	double THD = TMath::Sqrt( sqSumHarm ) / Arms;
	//		</THD>

	//		<ALTRO>
	double SINAD1 = Arms / NAD;
	double eta = TMath::Sqrt( NAD*NAD - Arms*Arms*THD*THD );
	double SNR = Arms / eta;
	double ENOB1 = nBit - TMath::Log2(   NAD / (  1 / TMath::Sqrt(12)  )   );		//bisognava mettere 1 invece di LSB

		//--------------------------------------------------	
		//My Test:
	double Ak1=0;
	double Ak2=0;
	double maxadc = TMath::MaxElement(gr1->GetN(),gr1->GetY());
	double minadc = TMath::MinElement(gr1->GetN(),gr1->GetY());
	//		for(int i=1;i<histMag2->GetNbinsX()/2.;++i) cout<<"bin "<<i<<": "<<histMag2->GetBinContent( i )<<endl;
	for ( int bin = 2; bin < maxBin-1; bin++ ) {	
	  Ak1 += histMag2->GetBinContent( bin ) * histMag2->GetBinContent( bin );
	}

	for ( int bin = maxBin+1; bin < histMag2->GetNbinsX()/2.; bin++ ) {	
	  Ak2 += histMag2->GetBinContent( bin ) * histMag2->GetBinContent( bin );
	}
	double Am=  histMag2->GetBinContent( maxBin );
	double mySinad=10*TMath::Log10((Am*Am)/(Ak1+Ak2));
	cout<<"max bin: "<<maxBin<<endl;
	cout<<"NBins: "<<histMag2->GetNbinsX()<<endl;
	cout<<"Am "<<Am<<" Ak1 "<<Ak1<<" Ak2 "<<Ak2<<endl;
	cout<<"mySinad: "<<mySinad<<endl;
	cout<<"myEnob "<< (mySinad - 1.76 + 20*TMath::Log(4095./((maxadc-minadc))))/6.02<<endl;

	//		</ALTRO>

	//	</ANALISI>

	//	<OUTPUT>

	//		<GRAFICO>

	if (grafica || salva) {



		histMagTot->GetXaxis()->SetRangeUser( 0, nPts/2 + 1 );
		histMagTot->GetXaxis()->SetTitle("Frequency bin");
		histMagTot->GetYaxis()->SetTitle("Amplitude");
		histMagTot->SetTitle("|X(#nu)|");
		histMagTot->SetLineColor(kRed + 2);

		histMag2->SetTitle("|X(#nu)|/|X_{in}|");
		histMag2->GetXaxis()->SetTitle("Frequency[Hz]");
		histMag2->GetYaxis()->SetTitle("Normalized amplitude");
		histMag2->SetLineColor(kRed + 2);
		histMag2->Scale( 1. / histMagTot->GetBinContent(maxBin) );
		//--------------------------------------------------	
		//My Test:
		double Ak1=0;
		double Ak2=0;
		double maxadc = TMath::MaxElement(gr1->GetN(),gr1->GetY());
		double minadc = TMath::MinElement(gr1->GetN(),gr1->GetY());
		//		for(int i=1;i<histMag2->GetNbinsX()/2.;++i) cout<<"bin "<<i<<": "<<histMag2->GetBinContent( i )<<endl;

		double Am=  histMag2->GetBinContent( maxBin );
		
		//--------------------------------------------------
		
		histMag2->GetXaxis()->SetLimits( 0, fSamp );
		histMag2->GetXaxis()->SetRangeUser(0, fSamp/2);
		

		if ( split == 0 )histMag3->SetTitle("FFT amplitude");
		if ( split == 1 )histMag3->SetTitle("FFT amplitude, odd samples");
		if ( split == 2 )histMag3->SetTitle("FFT amplitude, even samples");
		if ( split == 3 )histMag3->SetTitle("FFT amplitude, 4by4 samples");

		histMag3->GetXaxis()->SetTitle("Frequency[Hz]");
		histMag3->GetYaxis()->SetTitle("20 #upoint log_{10}(X(#nu)/X_{0}) [dB]");
		histMag3->SetLineColor(kRed + 2);
		histMag3->Scale( 1 / histMagTot->GetBinContent(maxBin) );
		histMag3->GetXaxis()->SetLimits( 0, fSamp );
		histMag3->GetXaxis()->SetRangeUser(0, fSamp/2);

		//plot logaritmico
		for (int bin = 1; bin <= histMag3->GetNbinsX(); bin++ ){

			if ( histMag2->GetBinContent(bin) == 0 ){

			  	histMag3->SetBinContent(bin, -2500);	//valore negativo arbitrario

			} else histMag3->SetBinContent(  bin, 20* TMath::Log10( histMag2->GetBinContent(bin) )  );

		}

	}

	if ( 1==1 ){

	  //		TCanvas * c1 = new TCanvas("c1", "", 1200, 400);
	  //		TCanvas * c2 = new TCanvas("c2", "", 1200, 400);
	  //	TCanvas * c3 = new TCanvas("c3", "", 1200, 400);
		TCanvas * c4 = new TCanvas("c4", "", 1200, 400);

		//		if (!minimalPlotting) c1->Divide(4, 1);

		if (!minimalPlotting){



//		c2->cd();
//		gStyle->SetOptStat(10);
//		gPad->SetGrid();
//		gPad->SetLogy();
//		histMagTot->DrawCopy("HIST");
//
//		c3->cd();
//		gStyle->SetOptStat(10);
//		gPad->SetGrid();
//		gPad->SetLogy();
//		histMag2->DrawCopy("hist");
//
			c4->cd();
			gStyle->SetOptStat(10);
			gPad->SetGrid();
			histMag3->SetLineWidth(2);
			histMag3->SetStats(0);
			histMag3->SetTitle("");
			histMag3->DrawCopy("hist");
			

		} else {

			gStyle->SetOptStat(10);
			gPad->SetGrid();
			histMag3->SetStats(0);
			histMag3->SetTitle("");
			histMag3->SetLineWidth(2);
			histMag3->DrawCopy("hist");

		}

	}

	//		</GRAFICO>

	//		<TESTUALE>

	if (output){

		std::cout << "============================================ Results ============================================" << std::endl;

		std::cout << "Dataset: " << nomeFile.c_str() << std::endl;
		std::cout << "Type: " << sampleType.c_str() << std::endl;
		std::cout << "Number of samples: " << nPts << std::endl;
		std::cout << "Input signal frequency [Hz]: " << histMag2->GetBinLowEdge(maxBin) << " +- " <<  histMag2->GetBinWidth(maxBin) << std::endl;

		std::cout << "NAD [ch] = " << NAD << std::endl;
		std::cout << "THD [%] = " << THD * 100 << std::endl;	//controllare che vada moltiplicato per 100, sempre che serva
		std::cout << "Arms [ch] = " << Arms << std::endl;
		std::cout << "SNDR [dB] = " << SINAD1 << " = " << 20 * TMath::Log10( SINAD1 ) << std::endl;
		//~ std::cout << "SINAD2 = " << SINAD2 << std::endl;
		//		std::cout << "SNR [dB] = " << SNR << " = " << 20  * TMath::Log10( SNR ) << std::endl;
		double a = 20  * TMath::Log10( SNR );
		cout<<"a "<<a<<endl;
		cout<<"freq "<<setprecision(10)<<freq<<endl;
		double jitter = TMath::Power(10,(-a/20.))/(2*3.14*freq);
		cout<<"!!!!!!!!!!!JITTER!!!!!!!"<< jitter<<endl;
		std::cout << "ENOB [bit] = " << ENOB1 << std::endl;
		//~ std::cout << "ENOB2 = " << ENOB2 << std::endl;

		std::cout << "=================================================================================================" << std::endl;

	}

	//		</TESTUALE>

	//		<SALVATAGGIO>

	if (salva){

		Risultati * res = (Risultati *) tf1.Get("Risultati");
		if ( res == nullptr ) res = new Risultati();
		res->AddFftData("Arms", Arms);
		res->AddFftData("NAD", NAD);
		res->AddFftData("THD", THD);
		res->AddFftData("ENOB", ENOB1);
		res->AddFftData( "SNDR", 20 * TMath::Log10(SINAD1) );
		res->AddFftData("SNR",  20 * TMath::Log10(SNR) );
		res->AddFftData( "freq", histMag2->GetBinCenter(maxBin-1) );
		res->AddFftData( "errFreq", histMag2->GetBinWidth(maxBin) ) / 2;
		res->Write("RisultatiFFT", TObject::kOverwrite);
		delete res;


		histMagTot->SetName("histMag");
		histMagTot->Write("", TObject::kOverwrite);
		histMag2->SetName("histMag2");
		histMag2->Write("", TObject::kOverwrite);
		histMag3->SetName("histMag3");
		histMag3->Write("", TObject::kOverwrite);
		
	}

	//		</SALVATAGGIO>

	// </OUTPUT>

	// <PULIZIA>

	delete histMagTot, histMag2;
	delete line1, line2;

	tf1.Close();

	// </PULIZIA>

	return 0;

}

	//

	//**sullo standard c'è scritto di prendere i bin n_i e ( nPts - n_i ), ma il contenuto è lo stesso. Controllo
	//~ std::cout << "Contenuto maxBin: " << histMag->GetBinContent( maxBin ) << std::endl;
	//~ std::cout << "Contenuto ( nPts - maxBin + 2 ): " << histMag->GetBinContent( nPts - maxBin + 2) << std::endl;
	// nPts - n_i = nPts - (maxBin - 1) = nPts - maxBin + 1 <--- questo è il bin di frequenza, non il bin dell'istogramma, per quello
	//devo ancora sommare 1: nPts - maxBin + 2

	//

	//	*NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
	//(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
	//La prima frase dice come passare da bin di frequenza a frequenza (proprio in Hz). Riscalare l'asse y mi serve
	//	=======> nello standard c'è una convenzione diversa: viene messo un fattore 1 / N alla DFT inversa, mentre la DFT non viene normalizzata

	//

	//~ double xMin = gPad->GetUxmin();
	//~ double xMax = gPad->GetUxmax();
	//~ double yMin = gPad->GetUymax() + 1500000;
	//~ double yMax = gPad->GetUymax() + 1500000;


	//~ std::cout << gPad->GetUxmax() << std::endl;
	//~ std::cout << gPad->GetUymax() << std::endl;

	//~ TGaxis * tg1 = new TGaxis( xMin, yMin, xMax, yMax, 0, fSamp/2 , 20 + 100*20, "-");
	//~ tg1->SetTitle("Frequenza[Hz]");
	//~ tg1->DrawClone();

	//~ grImFFT.DrawClone();

	//~ c1->cd(4);
	//~ gPad->SetGrid();
	//~ histPh->GetXaxis()->SetRangeUser( 0, nPts/2 + 1 );
	//~ histPh->GetXaxis()->SetLimits( 0, MSPS*1000000 );
	//~ histPh->GetXaxis()->SetLimits( 0, 1 );
	//~ histPh->GetXaxis()->SetTitle("Frequenza[Hz]");
	//~ histPh->GetXaxis()->SetTitle("Bin di frequenza");
	//~ histPh->GetYaxis()->SetTitle("Fase");
	//~ histPh->SetTitle("Fase della DFT");
	//~ histPh->SetLineColor(kRed + 2);
	//~ histPh->DrawCopy();
	//~ grReFFT.DrawClone();

	//~ c1->cd(4);
	//~ gStyle->SetOptStat(10);
	//~ gPad->SetGrid();
	//~ gPad->SetLogy();
	//~ histMag->GetXaxis()->SetTitle("Frequenza[Hz]");
	//~ histMag->GetXaxis()->SetLimits( 0, fSamp );
	//~ histMag->DrawCopy();

	//
