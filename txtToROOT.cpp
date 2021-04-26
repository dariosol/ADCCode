
//~ #include "Header.h"

#if !defined (__CLING__) || defined (__ROOTCLING__)
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <math.h>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TPad.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>		//per usare atof
#include <vector>
#include <bitset>
#include <chrono>

#endif

void txtToROOT(

	       /*	std::string nomeFile = "Dati/dati.txt",*/
	       std::string nomeFile = "~/LiteDTU/Files/02_08_board2_adctm_sin500k_clk160_800400-2.csv.txt",
	bool save = 0,
	bool print = 1,
	bool binaryMode = 1,
	bool split = 0		//divide i campionamenti pari e dispari

	){
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	//definizioni varie
	std::string nomeIst = "histSamp";
	std::string nomeIstEven = nomeIst + "Even";
	std::string nomeIstOdd = nomeIst + "Odd";
	std::string nomeIst4by4 = nomeIst + "FourByFour";
	
	std::string nomeTGraph = "grSamp";
	std::string nomeTGraphEven = nomeTGraph + "Even";
	std::string nomeTGraphOdd = nomeTGraph + "Odd";
	std::string nomeTGraph4by4 = nomeTGraph + "FourByFour";
	
	std::string nomeOutput = nomeFile;
	for (int i = 0; i < 4; i++) nomeOutput.pop_back();		//assumendo ".txt"
	nomeOutput.append(".root");

	//definizioni iniziali parametri ADC
	int nBit = 12;
	int nCh = (int)TMath::Power(2, nBit);
	int chMax = nCh - 1;

	//istogramma frequenze
	TH1D hist( nomeIst.c_str(), "Distribution of ADC codes", nCh, -0.5, chMax + 0.5 );
	hist.GetXaxis()->SetTitle("ADC code");
	hist.GetYaxis()->SetTitle("Frequency");
	hist.SetFillColor(kRed + 2);
	hist.SetLineColor(kRed + 2);

	//istogramma frequenze, campionamenti pari
	TH1D histEven( nomeIstEven.c_str(), "Distribution of ADC codes, even samples", nCh, -0.5, chMax + 0.5 );
	histEven.GetXaxis()->SetTitle("ADC code");
	histEven.GetYaxis()->SetTitle("Frequency");
	histEven.SetFillColor(kRed + 2);
	histEven.SetLineColor(kRed + 2);

	//istogramma frequenze, campionamenti dispari
	TH1D histOdd( nomeIstOdd.c_str(), "Distribution of ADC codes, odd samples", nCh, -0.5, chMax + 0.5 );
	histOdd.GetXaxis()->SetTitle("ADC code");
	histOdd.GetYaxis()->SetTitle("Frequency");
	histOdd.SetFillColor(kRed + 2);
	histOdd.SetLineColor(kRed + 2);

	//istogramma frequenze, campionamenti dispari
	TH1D hist4by4( nomeIst4by4.c_str(), "Distribution of ADC codes, odd samples", nCh, -0.5, chMax + 0.5 );
	hist4by4.GetXaxis()->SetTitle("ADC code");
	hist4by4.GetYaxis()->SetTitle("Frequency");
	hist4by4.SetFillColor(kRed + 2);
	hist4by4.SetLineColor(kRed + 2);

	//apro il file di testo
	std::ifstream fileDati;
	fileDati.open( nomeFile.c_str() );

	std::string riga = "";

	long i = 0;
	long iEven = 0;
	long iOdd = 0;
	long i4by4 = 0;
	
	int adc = 0;

	TGraph tGr;
	tGr.SetName( nomeTGraph.c_str() );
	tGr.SetTitle("ADC codes vs sample number");
	tGr.SetMarkerStyle(25);
	tGr.GetXaxis()->SetTitle("Sample number");
	tGr.GetYaxis()->SetTitle("ADC code");

	TGraph tGrEven;
	tGrEven.SetName( nomeTGraphEven.c_str() );
	tGrEven.SetTitle("ADC codes vs even sample number");
	tGrEven.SetMarkerStyle(25);
	tGrEven.GetXaxis()->SetTitle("Sample number");
	tGrEven.GetYaxis()->SetTitle("ADC code");

	TGraph tGrOdd;
	tGrOdd.SetName( nomeTGraphOdd.c_str() );
	tGrOdd.SetTitle("ADC codes vs odd sample number");
	tGrOdd.SetMarkerStyle(25);
	tGrOdd.GetXaxis()->SetTitle("Sample number");
	tGrOdd.GetYaxis()->SetTitle("ADC code");

	TGraph tGr4by4;
	tGr4by4.SetName( nomeTGraph4by4.c_str() );
	tGr4by4.SetTitle("ADC codes vs odd sample number");
	tGr4by4.SetMarkerStyle(25);
	tGr4by4.GetXaxis()->SetTitle("Sample number");
	tGr4by4.GetYaxis()->SetTitle("ADC code");
	
	Header header;
	int sample=0;
	
	//ciclo sul file di testo
	if ( fileDati.is_open() ){
	  std::cout<<"file opened "<<std::endl;
		while ( std::getline(fileDati, riga) ){

			int stringSize = riga.size();
			
			if ( riga.find("#") < stringSize ) {
			  
				//estrazione key e value
				while ( riga.find(" ") < stringSize ) riga.erase( riga.find(" ") , 1 );	//rimozione degli spazi
				riga.erase( riga.find("#") , 1 );

				int pos = riga.find("=");
				std::string key = riga.substr(0, pos);
				std::string valString = riga.substr( pos + 1, stringSize - pos );
				double value = std::stod(valString);
				std::cout<<"key in the Header: "<<key<<" value "<<value<<std::endl;
				header.AddInfo( key, value );

				continue;

			}

			double adcRaw = -1;

			try{

				adcRaw = std::stod( riga );

			} catch ( const std::invalid_argument& invArg ) {

				std::cerr << adcRaw << " ---> " << "Argomento invalido a " << invArg.what() << " nella riga " << i + 1 << std::endl;
				continue;

			}

			if (!binaryMode) {

				adc = round( adcRaw );

			} else {

				adc = (int) std::strtoul( riga.c_str(), nullptr, 2 );

			}
			
			//cout<<"sample "<<i<<" ADC: "<<adc<<endl;
			sample++;
			
			tGr.SetPoint( i, i, adc );

			hist.Fill( adc );

			if ( i % 2 == 0 ){

				tGrEven.SetPoint( iEven, iEven, adc );
				histEven.Fill( adc );

				iEven++;

			} else {

				tGrOdd.SetPoint( iOdd, iOdd, adc );
				histOdd.Fill( adc );

				iOdd++;

			}

			if(i % 4 ==0) {
			  tGr4by4.SetPoint( i4by4, i4by4, adc );
			  hist4by4.Fill( adc );
			  i4by4++;
			}

			i++;	//<==== controllare questo

		}

	}



	TH1D * ADC_ph0 = new TH1D( "ADC_ph0", "adc_ph0", 1000, -0.5, 999.5 );
	TH1D * ADC_ph1 = new TH1D( "ADC_ph1", "adc_ph1", 1000, -0.5, 999.5 );
	TH1D * ADC_ph2 = new TH1D( "ADC_ph2", "adc_ph2", 1000, -0.5, 999.5 );
	
	bool dohistoPhase=1;
	int cycle=-1;
	int j=0;
	
	hist.DrawCopy();
	cout<<"Graph points: "<< tGr.GetN()<<endl;

       	Double_t *y  =tGr.GetY();
	auto dataSize =  tGr.GetN();
	cout<<"!!!size!!! "<<dataSize<<endl;

	vector<int> newy;
	vector<int> newx;
	int k=0;
	int oldvalue=0;
	TH1F * DeltaADC = new TH1F("DeltaADC","",200,-100,100);
	TH1F * DeltaADC_New = new TH1F("DeltaADC_New","",200,-100,100);
	int removed=0;
	for(int i=1; i<dataSize;++i) {
	  DeltaADC->Fill(y[i]-y[i-1]);
	  if(fabs(y[i]-oldvalue) < 60) {
	    if(removed!=0)cout<<"Removed "<<removed<<endl;
	    newy.push_back(y[i]);
	    newx.push_back(k);
	    k++;
	    oldvalue=y[i];
	    removed=0;
	  }
	  else {
	    removed++;
	    // cout<<"Removing sample "<<i<< " with value "<< y[i]<<endl;
	  }
	   
	}
	int* newx_arr = newx.data();
	int* newy_arr = newy.data();

	TGraph *tGrNew = new TGraph(newx.size(),newx_arr,newy_arr);
	tGrNew->SetName( "New" );
	tGrNew->SetTitle("ADC codes vs sample number");
	tGrNew->SetMarkerStyle(25);
	tGrNew->GetXaxis()->SetTitle("Sample number");
	tGrNew->GetYaxis()->SetTitle("ADC code");
	Double_t *ynew  =tGrNew->GetY();

	TH1D histNew( "histNew", "Distribution of ADC codes New", nCh, -0.5, chMax + 0.5 );
	histNew.GetXaxis()->SetTitle("ADC code new");
	histNew.GetYaxis()->SetTitle("Frequency");
	histNew.SetFillColor(kRed + 2);
	histNew.SetLineColor(kRed + 2);

	for(int i=0; i< tGrNew->GetN();++i) {
	  histNew.Fill( ynew[i] );
	  DeltaADC_New->Fill(ynew[i]-ynew[i-1]);
	  if((i%262000)==0) {
	    dohistoPhase=1;
	    cycle++;
	  }
	  if(dohistoPhase) {
	    cout<<j<<" "<< ynew[i]<<endl;
	    if(cycle==0) ADC_ph0->SetBinContent( j+1, ynew[i]);
	    if(cycle==1) ADC_ph1->SetBinContent( j+1, ynew[i]);
	    if(cycle==2) ADC_ph2->SetBinContent( j+1, ynew[i]);
	    j++;
	    if(j==1000) {
	      dohistoPhase=0;
	      j=0;
	    }
	  }

	}





	
	const int nPts = i;
	const int nEven = iEven;
	const int nOdd = iOdd;
	const int n4by4 = i4by4;
	
	fileDati.close();

	if (1==1) {

		//rette per campionamenti
		TF1	* line1 = new TF1( "line1", "0", 0, nPts - 1 );
		line1->SetLineColor(kOrange - 3);
		line1->SetLineWidth(1);
		line1->SetLineStyle(1);
		TF1	* line2 = new TF1( "line2", std::to_string(nCh).c_str(), 0, nPts - 1 );
		line2->SetLineColor(kOrange - 3);
		line2->SetLineWidth(1);
		line2->SetLineStyle(1);

		//grafica
		TCanvas * c1 = new TCanvas("c1");
		//	c1->Divide(2,1);
		//	c1->GetPad(1)->SetLeftMargin(0.15);			//serve per far sì che le label degli assi Y non
		//	c1->GetPad(2)->SetLeftMargin(0.15);			//vengano tagliate perché il margine è troppo piccolo

		//	c1->cd(1);
		gPad->SetGrid();
		hist.DrawCopy();

		//	c1->cd(2);
		TCanvas * c2 = new TCanvas("c2");

		gPad->SetGrid();
		tGr.SetMarkerStyle(20);
		tGr.SetMarkerSize(0.5);
		tGr.GetYaxis()->SetRangeUser(-100, nCh + 100);
		tGr.DrawClone("AP");
		line1->DrawCopy("same");
		line2->DrawCopy("same");

		delete line1;
		delete line2;

	}

	if (save){

		//salvataggio
		TFile tf1( nomeOutput.c_str(), "recreate");
		hist.Write();
		tGr.Write();
		tGrNew->Write();
		histNew.Write();
		DeltaADC->Write();
		DeltaADC_New->Write();
		
		ADC_ph0->Write("", TObject::kOverwrite);
		ADC_ph1->Write("", TObject::kOverwrite);
		ADC_ph2->Write("", TObject::kOverwrite);
		if ( split ){

			tGrEven.Write();
			tGrOdd.Write();
			tGr4by4.Write();

			histEven.Write();
			histOdd.Write();
			hist4by4.Write();

		}

		header.Write();
		tf1.Close();

	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "[txtToROOT]:Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
 
}
