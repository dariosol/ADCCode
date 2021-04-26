
#include "Risultati.h"
#include "Header.h"

#if !defined (__CLING__) || defined (__ROOTCLING__)
//		C++
#include <cmath>
#include <iostream>
#include <unordered_map>
//		ROOT
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>
#include <TStyle.h>
#include <TH2D.h>
#endif

//NOTA BENE: se non si conosce bene l'ampiezza dell'onda allora l'SNDR è inaccurato!

//	<DA FARSI>
//	</DA FARSI>

int fit(

	std::string nomeFile = "datiFFiT.root",

	bool salva = 0,
	Double_t MSPS = 160,
	bool useCalibratedSamples = 1,		//analizza i dati calibrati, se presenti
	string ReportName = "ReportBoard5_16_09_2019.txt",
	bool report =0,
	int split = 0 		//0: tutti, 1: dispari, 2: pari

	){

	//definizioni iniziali parametri ADC
	const Double_t nBit = 12;
	const int nCh = (int) TMath::Power(2, nBit);
	const int chMax = nCh - 1;
	const Double_t vmax = 0.6;
	const Double_t vmin = -0.6;
	const Double_t FSR = vmax - vmin;
	const Double_t LSB = FSR / nCh;
	//~ const Double_t MSPS = 160;
	if ( split ==1 or split ==2 ) MSPS = MSPS / 2.;
	if ( split ==3 ) MSPS = MSPS / 4.;
	
	const Double_t fSamp = MSPS * 1000000;
	const Double_t tSamp = 1 / fSamp;

	//definizioni varie
	//~ std::string nomeTGraph = "tGr";
	const std::string nomeTGraphDefault = "grSamp";
	std::string nomeTGraph = nomeTGraphDefault;
	if ( split == 1 && useCalibratedSamples == 0 ) nomeTGraph += "Odd";
	if ( split == 2 && useCalibratedSamples == 0 ) nomeTGraph += "Even";
	if ( split == 3 && useCalibratedSamples == 0 ) nomeTGraph += "FourByFour";

	if ( useCalibratedSamples ) nomeTGraph += "Calib";
	const std::string ampKey = "amp[V]";
	const std::string freqKey = "freq[Hz]";
	std::string sampleType = "";
	if ( !useCalibratedSamples && !split ) sampleType += "All";
	if ( split == 1 ) sampleType += "Odd";
	if ( split == 2 ) sampleType += "Even";
	if ( split == 3 ) sampleType += "FourByFour";
	if ( useCalibratedSamples ) sampleType += "Calibrated";

	const bool grafica = 1;
	const bool output = 1;

	//	<LETTURA FILE>

	TFile tf1( nomeFile.c_str(), "update");

	//controllo esistenza del file
	if ( !tf1.IsOpen() ){

		std::cout << "File \"" << nomeFile.c_str() << "\" not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		return 1;

	}

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
		return 2;

	}

	gr1->SetMarkerStyle(7);

	//		<LETTURA HEADER>
	Header * headerObj = (Header *) tf1.Get( "Header" );

	if ( headerObj == nullptr ){

		std::cout << "Header not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete gr1;
		tf1.Close();
		return 3;

	}

	std::unordered_map < std::string, Double_t > headerMap = headerObj->GetHeader();
	delete headerObj;

	Double_t amp = 0;
	Double_t freq = 0;

	try {

		amp = headerMap.at( ampKey );
		freq = headerMap.at( freqKey );
		std::cout<<"amp: "<<amp<<std::endl;
		std::cout<<"freq: "<<freq<<std::endl;

	} catch ( const std::out_of_range& invArg ) {

		std::cout << "Key not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete gr1;
		tf1.Close();
		return 4;

	}

	//		</LETTURA HEADER>
	
	//altre definizioni
	Double_t omega = 2 * TMath::Pi() * freq;
	Double_t T = 1 / freq;
	Double_t Tsamp = T / tSamp;	//periodo in unità di campionamenti ideali
	Double_t omegaSamp = 2 * TMath::Pi() / Tsamp;

	int nPts = gr1->GetN();
	//int nPts = 16384;

	Double_t * iSamp = new Double_t[ nPts ];
	Double_t * res = new Double_t[ nPts];
	Double_t * ch = new Double_t[ nPts ];

	//	</LETTURA FILE>
	cout<<"Frequency set by user: "<<setprecision(10)<<freq<<endl;
	//rette
	TF1	* line1 = new TF1( "line1", "0.5", 0, nPts - 1 );
	line1->SetLineColor(kOrange - 3);
	line1->SetLineWidth(1);
	line1->SetLineStyle(1);

	TF1	* line2 = new TF1( "line2", "-0.5", 0, nPts - 1 );
	line2->SetLineColor(kOrange - 3);
	line2->SetLineWidth(1);
	line2->SetLineStyle(1);

	TF1	* line3 = new TF1( "line3", "0", 0, nPts - 1 );
	line3->SetLineColor(kOrange - 3);
	line3->SetLineWidth(1);
	line3->SetLineStyle(1);

	TF1	* line4 = new TF1( "line4", std::to_string(nCh).c_str(), 0, nPts - 1 );
	line4->SetLineColor(kOrange - 3);
	line4->SetLineWidth(1);
	line4->SetLineStyle(1);

	//istogramma distribuzione residui
	TH1F histRes( "histRes", "histRes", 800, -20, 20);

	//istogramma campionamenti in funzione della fase
	TH1D histPhase( "histPhase", "Sampling phase", nPts, 0 - 2 * TMath::Pi() / ( 2 * nPts ) , 2 * TMath::Pi() + 2 * TMath::Pi() / ( 2 * nPts ));

	//scatter plot dei residui in funzione della fase
	TH2D histErrPhase("histErrPhase", "Phase-residuals distribution", 50, 0, 2*TMath::Pi(), 5000, -50, 50 );

	//ricava valori iniziali dei parametri del fit dai dati
	int nMin_exp = TMath::LocMin( nPts, gr1->GetY() );
	int nMax_exp = TMath::LocMax( nPts, gr1->GetY() );
	Double_t chMin_exp = gr1->GetY()[ nMin_exp ];
	Double_t chMax_exp = gr1->GetY()[ nMax_exp ];
	Double_t A_exp = ( chMax_exp - chMin_exp ) / 2;
	Double_t offset_exp = chMax_exp - A_exp;

//fit
	TF1 * sinFit = new TF1( "sinFit", "[0] * TMath::Sin( [1] * x + [2] ) + [3]", 0, nPts -1  );
	sinFit->SetParameter(0, A_exp);
	sinFit->SetParLimits(0, A_exp*0.95, A_exp*1.05);//Prima 0.7, 1.3
	sinFit->SetParName(0, "Amplitude");
	//sinFit->FixParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	sinFit->SetParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	sinFit->SetParLimits(1, omegaSamp * 0.99, omegaSamp * 1.01);
	sinFit->SetParName(1, "Omega");
	//~ sinFit->SetParameter(2, phi);
	sinFit->SetParName(2, "Phase");
	sinFit->SetParameter(3, offset_exp);
	sinFit->SetParLimits(3, offset_exp*0.95, offset_exp*1.05);
	sinFit->SetParName(3, "Offset");
	//~ gr1->Fit("sinFit", "RN0Q");		//fithttps://stackoverflow.com/questions/160930/how-do-i-check-if-an-integer-is-even-or-odd
	gr1->Fit("sinFit", "RN0QB");		//fit

	//calcolo varianza
	Double_t var = 0;
	int broken=0;
	//	TH1D *histINL =new TH1D("histINL", "INL", 4096, -0.5, 4096 - 0.5);
	//	TH1D *histNorm=new TH1D("histNorm", "", 4096,-0.5, 4096 - 0.5);
	//TH1D *histMissingCodes=new TH1D("histMissingCosdes", "", 4096,-0.5, 4096 - 0.5);

	sinFit->Print();
	for ( long i = 0; i < nPts; i++ ) {

		Double_t x, y = 0;
		gr1->GetPoint(i, x, y);
		iSamp[i] = x;	
		res[i] = y - sinFit->Eval( x );		
		var += res[i] * res[i];
		histRes.Fill( res[i] );
		//		histINL->AddBinContent( y+1, res[i] );

		//		histNorm->AddBinContent( y+1, 1 );	
		//fase
		Double_t fase = fmod( sinFit->GetParameter(1) * i + sinFit->GetParameter(2),  2 * TMath::Pi() );
		histPhase.Fill( fase );	
		//scatter plot
		histErrPhase.Fill(fase, res[i]);
	}

	//	cout<<"********before "<<histINL->GetBinContent(2500)<<endl;
	
//	histINL->Divide(histNorm);
//	for(int i=0;i<4096;++i) {
//	  if(histINL->GetBinContent(i)!=0) histINL->SetBinError(i,histINL->GetBinContent(i)/TMath::Sqrt(histNorm->GetBinContent(i)));
//	}
//	
//
//	TH1D *histDNL =new TH1D("histDNL", "DNL", 4096, -0.5, 4096 - 0.5);
//
//	for(int i=1; i < 4096; i++) {
//	  if(histINL->GetBinContent(i+1)==0) continue;
//	  histDNL->SetBinContent(i,histINL->GetBinContent(i+1) - histINL->GetBinContent(i));
//	  histDNL->SetBinError(i,0.1);
//	}
//
//	for(int i=2; i < 4096; i++) {
//	  if(histNorm->GetBinContent(i+1)==0) {
//	    histMissingCodes->Fill(i);
//	    }
//	}
//	
	//grafico dei residui in funzione del tempo
	TGraph gr2( nPts, iSamp, res );
	gr2.SetName( "grResT" );
	TF1 *fg = new TF1("fg","gaus",-10,10);
	cout<<"*********"<<endl;
	histRes.Fit(fg);
	
	cout<<"*********"<<endl;
	cout<<"par 0 "<<fg->GetParameter(0)<<endl;
	cout<<"par 1 "<<fg->GetParameter(1)<<endl;
	cout<<"par 2 "<<fg->GetParameter(2)<<endl;
	
	//grafico dei residui tra fit e valore vero
	// TGraph gr3( nPts, t, res2);
	// gr3.SetMarkerStyle(6);
	// gr3.SetLineStyle(3);

	//varianze 
	var = var * LSB * LSB / (nPts);
	Double_t varQ = ( LSB * LSB ) / 12;
	Double_t sigmaQ = TMath::Sqrt( varQ );
	Double_t NAD = TMath::Sqrt(var);
	//ENOB:
	Double_t enob1 = nBit - TMath::Log2( NAD / sigmaQ );	
	Double_t enob3 = TMath::Log2(  FSR / ( TMath::Sqrt(12) * (fg->GetParameter(2))*(LSB)));

	Double_t errENOB =(1./(fg->GetParameter(2)*TMath::Log(2)))*fg->GetParError(2);
	cout<<"NAD: "<<NAD<<endl;
	cout<<"sigma: "<<(fg->GetParameter(2)*LSB)<<endl;
	cout<<"*******ENOB 1 "<<enob1<<endl;
	cout<<"*******ENOB 3 "<<enob3<<" pm "<<errENOB<<endl;
	//~ Double_t enob3 = nBit - TMath::Log2( histRes.GetRMS()*LSB / sigmaQ );
	
	//SNDR:
	//Double_t Arms = amp / TMath::Sqrt(2); //amp e' un valore a caso
	Double_t Arms=sinFit->GetParameter(0)*LSB / TMath::Sqrt(2);
	Double_t sinad1 = 20 * TMath::Log10( Arms / NAD );
	Double_t maxadc = TMath::MaxElement(gr1->GetN(),gr1->GetY());
	Double_t minadc = TMath::MinElement(gr1->GetN(),gr1->GetY());

	std::cout<<"max adc: "<<maxadc<<" min adc: "<<minadc<<std::endl;
	Double_t sinad2 = ((enob1 * 6.02) + 1.76 - 20*TMath::Log(4095./(sinFit->GetParameter(0)+sinFit->GetParameter(3)))); //if full scale aplitude.
	cout<<"sinad1: "<<sinad1<<endl;
	cout<<"sinad2: "<<sinad2<<endl;
	
	if (output){
	  ofstream outfile;
	  if(report)outfile.open(ReportName, std::ofstream::app);
		std::cout << "======================================= Results =======================================" << std::endl;

		std::cout << "Dataset: " << nomeFile.c_str() << std::endl;
		if(report) outfile << "Dataset: " << nomeFile.c_str() << std::endl;
		
		std::cout << "Type: " << sampleType.c_str() << std::endl;
		if(report) outfile << "Type: " << sampleType.c_str() << std::endl;
		
	       	std::cout << "Effective Resolution: " << TMath::Log2(FSR/NAD) << std::endl;
	       	if(report) outfile << "Effective Resolution: " << TMath::Log2(FSR/NAD) << std::endl;
	       
		std::cout << "Sine fit parameters:" <<  sinFit->GetParameter(0)<<"*sin("<< sinFit->GetParameter(1)<<"*x + "<< sinFit->GetParameter(2)<<") + "<<sinFit->GetParameter(3)<<std::endl;
	       	if(report) outfile << "Sine fit parameters:" << std::endl;
		
		std::cout << "\t0) Amplitude: [ch]" << sinFit->GetParameter(0) << std::endl;
		if(report) outfile << "\t0) Amplitude: " << sinFit->GetParameter(0) << std::endl;
		
		std::cout << "\t1) Frequency [Hz]: " <<setprecision(10)<<sinFit->GetParameter(1) / ( 2 * TMath::Pi() ) * fSamp << std::endl;	//è in unità di campionamenti
		if(report) outfile << "\t1) Frequency [Hz]: " << sinFit->GetParameter(1) / ( 2 * TMath::Pi() ) * fSamp << std::endl;	//è in unità di campionamenti

		std::cout << "\t2) Phase: " << sinFit->GetParameter(2) << std::endl;
		if(report) outfile << "\t2) Phase: " << sinFit->GetParameter(2) << std::endl;

		std::cout << "\t3) Offset: " << sinFit->GetParameter(3) << std::endl;
		if(report) outfile << "\t3) Offset: " << sinFit->GetParameter(3) << std::endl;

		//V in from amplitude:
		std::cout << "\t0) Amplitude+ [V]: " << sinFit->GetParameter(0)*(1.2/4095) + 0.6 << std::endl;
		std::cout << "\t0) Amplitude- [V]: " << sinFit->GetParameter(0)*(1.2/4095) - 0.6 << std::endl;
		std::cout << "\t0) Max channel expected: " << sinFit->GetParameter(0) + 2048 << std::endl;
		
		
		std::cout << "Variance [LSB]: " << var / ( LSB * LSB ) << std::endl;
		if(report) outfile << "Variance [LSB]: " << var / ( LSB * LSB ) << std::endl;

		std::cout << "Variance (ideal) [LSB]: " << varQ / ( LSB * LSB ) << std::endl;
		if(report) outfile << "Variance (ideal) [LSB]: " << varQ / ( LSB * LSB ) << std::endl;

		//~ std::cout << "ENOB, metodo 1: " << enob1 << std::endl;
		std::cout << "ENOB [bit]: " << enob1 << std::endl;
		if(report) outfile << "ENOB [bit]: " << enob1 << std::endl;
		
		//~ // std::cout << "ENOB, metodo 3: " << enob3 << std::endl;

		std::cout << "SNDR [dB]: " << sinad1 << std::endl;	//questo è il metodo presente nello standard IEEE 2017
		if(report) outfile << "SNDR [dB]: " << sinad1 << std::endl;	//questo è il metodo presente nello standard IEEE 2017

		std::cout << "SNDR - from ENOB, assuming full range input [dB]: " << sinad2 << std::endl;	//from ENOB
		if(report) outfile << "SNDR - from ENOB, assuming full range input [dB]: " << sinad2 << std::endl;	//from ENOB
		if(report) outfile<<std::endl;
		//~ std::cout << "SINAD[dB], direttamente da ENOB:" << sinad2 << std::endl;

		std::cout << "========================================================================================" << std::endl;

	}

	// <OUTPUT>
	
	if (grafica || salva){

		if (split == 0) gr1->SetTitle("Sampled signal");
		if (split == 1) gr1->SetTitle("Sampled signal, odd samples");
		if (split == 2) gr1->SetTitle("Sampled signal, even samples");
		if (split == 3) gr1->SetTitle("Sampled signal, 4by4 samples");
		
		if (split == 0) gr1->GetXaxis()->SetTitle("Sample");
		if (split == 1) gr1->GetXaxis()->SetTitle("Sample (odd)");
		if (split == 2) gr1->GetXaxis()->SetTitle("Sample (even)");
		if (split == 3) gr1->GetXaxis()->SetTitle("Sample (4by4)");

		gr1->GetYaxis()->SetTitle("ADC code");
		gr1->GetXaxis()->SetRangeUser( 0, nPts * 0.05);
		gr1->GetYaxis()->SetRangeUser( -100, nCh + 100);

		gr2.SetTitle("Residuals as a function of time");
		if (split == 0) gr2.GetXaxis()->SetTitle("Sample");
		if (split == 1) gr2.GetXaxis()->SetTitle("Sample (odd)");
		if (split == 2) gr2.GetXaxis()->SetTitle("Sample (even)");
		if (split == 3) gr2.GetXaxis()->SetTitle("Sample (4by4)");
		
		gr2.GetYaxis()->SetTitle("#Delta[ADC counts]");
		gr2.GetXaxis()->SetRangeUser( -1, nPts );
		gr2.SetLineStyle(1);
		gr2.SetLineWidth(1);
		gr2.SetMarkerStyle(6);
		// gr2.SetMarkerStyle(7);
		// gr2.SetMarkerSize(0.5);

		histRes.SetTitle( "Distribution of residuals" );
		histRes.GetXaxis()->SetTitle( "#Delta[ADC counts]" );
		histRes.GetYaxis()->SetTitle( "Counts" );
		histRes.GetXaxis()->SetRangeUser( histRes.GetMean() - 5*histRes.GetStdDev(), histRes.GetMean() + 5*histRes.GetStdDev() );
		//~ histRes.SetLineColor( kRed + 2 );
		//~ histRes.SetFillColor( kOrange - 3 );
		//~ histRes.SetFillStyle( 3004 );
		//~ histRes.GetXaxis()->SetNdivisions(20);

		histPhase.SetLineColor( kRed + 2 );
		//~ histPhase.SetFillColor( kOrange - 3 );
		//~ histPhase.SetFillStyle( 3004 );
		histPhase.GetXaxis()->SetTitle("Phase[rad]");
		histPhase.GetYaxis()->SetTitle("Counts");
		histPhase.SetMarkerStyle(7);

		histErrPhase.GetXaxis()->SetTitle("Phase[rad]");
		histErrPhase.GetYaxis()->SetTitle("#Delta[ADC counts]");
		histErrPhase.GetYaxis()->SetRangeUser( histRes.GetMean() - 5*histRes.GetStdDev(), histRes.GetMean() + 5*histRes.GetStdDev() );

		sinFit->SetLineWidth(1);
		sinFit->SetLineColor( kRed );

	}

	//if (grafica){
	  if (1==1) {
		TCanvas * c1 = new TCanvas("c1", "", 1250, 650);
		//	c1->Divide(2, 3);

		//grafico campionamenti
		//	c1->cd(1);
		gPad->SetGrid();
		gr1->SetMarkerStyle(20);
		gr1->DrawClone("AP");
		line3->DrawCopy("same");
		line4->DrawCopy("same");

		//grafico residui
		//		c1->cd(3);
		TCanvas * c2 = new TCanvas("c2", "", 1250, 650);

		gPad->SetGrid();
		gr2.DrawClone("LPA");
		line1->DrawCopy("same");
		line2->DrawCopy("same");

		//distribuzione dei residui
		TCanvas * c3 = new TCanvas("c3", "", 1250, 650);

		//c1->cd(4);
		gPad->SetGrid();
		histRes.DrawCopy();

		c3->Update();

		TCanvas * c4 = new TCanvas("c4", "", 1250, 650);
		//	c1->cd(5);
		gPad->SetGrid();
		gStyle->SetOptStat(0);
		histPhase.DrawCopy("");

		c4->Update();

		//scatter plot
		TCanvas * c5 = new TCanvas("c5", "", 1250, 650);

		//	c1->cd(6);
		gStyle->SetOptStat(0);
		gPad->SetGrid();
		histErrPhase.DrawCopy("COLZ");

		c5->Update();

		//grafico fit; deve essere l'ultimo!
		//	gr1->GetXaxis()->SetRangeUser(32680,32690);
		gr1->GetXaxis()->SetRangeUser(0,20);
		gr1->SetTitle("Sinewave fit");
		//		c1->cd(2);
		TCanvas * clast = new TCanvas("clast", "", 1250, 650);

		gPad->SetGrid();
		gStyle->SetOptStat(1110);
		gStyle->SetOptFit(1110);
		gr1->DrawClone("PA");
		sinFit->DrawCopy("same");
		line3->DrawCopy("same");
		line4->DrawCopy("same");
		gPad->Modified();

		//TCanvas * c6 = new TCanvas("c6", "", 1250, 650);
		//gPad->SetGrid();
		//gStyle->SetOptStat(0);
		//histINL->GetYaxis()->SetRangeUser(-4,4);
		//histINL->DrawCopy("e1p");
		//TLine *INLmin = new TLine(0,-1.5,4096,-1.5);
		//TLine *INLmax = new TLine(0,+1.5,4096,+1.5);
		//INLmin->SetLineColor(kRed);
		//INLmax->SetLineColor(kRed);
		//
		//INLmin->Draw("same");
		//INLmax->Draw("same");
		//
		//c6->Update();
		//
		//TCanvas * c7 = new TCanvas("c7", "", 1250, 650);
		//gPad->SetGrid();
		//gStyle->SetOptStat(0);
		//histDNL->SetMarkerSize(0.5);		
		//histDNL->SetMarkerStyle(20);
		//		histDNL->DrawCopy();
		//histMissingCodes->DrawCopy("");
		//c7->Update();
	}

	// <SALVATAGGIO>

	if (salva){

		Risultati * results = (Risultati *) tf1.Get("Risultati");
		if ( results == nullptr ) results = new Risultati();
		results->AddFitData( "Variance", var / ( LSB * LSB ) );
		results->AddFitData( "ENOB", enob1 );
		results->AddFitData( "SNDR", sinad1 );
		results->Write("Risultati", TObject::kOverwrite);
		delete results;

		gr2.Write("", TObject::kOverwrite);
		sinFit->Write("", TObject::kOverwrite);
		histRes.Write("", TObject::kOverwrite);
		histPhase.Write("", TObject::kOverwrite);
		histErrPhase.Write("", TObject::kOverwrite);

	}

	// </SALVATAGGIO>

	// </OUTPUT>

	delete sinFit;
	delete line1, line2, line3, line4;
	delete gr1;
	delete [] iSamp;
	delete [] ch;
	delete [] res;

	tf1.Close();

	return 0;

}
