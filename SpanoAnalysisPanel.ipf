#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

#include <KBColorizeTraces>

//imports a time series of UV-vis spectra acquired via the Ocean Optics QE-Pro spectrometer/software 
Function ITS(name,tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale, [d,load,doFits,wli,wlf,aggregateType])	//IIMPORTANT:Wavelengths must be in row 18 for proper data 
	String name //name of data set
	Variable tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale
	Variable d //display the imported data
	Variable load	//Load the data? If not, then the timewave,2D wave and lambda wave must aalready be present
	Variable doFits	//Run Spano analysis
	Variable wli,wlf
	String aggregateType
	
	if(ParamIsDefault(d))
		d = 1
	endif
	
	if(ParamIsDefault(load))
		load = 1
	endif
	
	if(ParamIsDefault(doFits))
		doFits = 1
	endif
	
	if(ParamIsDefault(wli))
		wli = 300
	endif
	
	if(ParamIsDefault(wlf))
		wlf = 800
	endif
	
	if(ParamIsDefault(aggregateType))
		aggregateType = "H"
	endif
	
	String CurrentFolder=GetDataFolder(1)
	String fname="root:"+name
	String wname
	NewDataFolder/o/s $fname
	
	//Load 2D data
	wname="C=2;N="+name+";"
	if(load)
		LoadWave/Q/J/M/U={1,2,1,2}/B=wname/K=0/L={12,13,0,1,0}/O/N/D
		if(V_Flag==0)
			Abort "No file selected"
		endif
		WAVE spec=$(name), lambda=$("CP_"+name), timeData=$("RP_"+name)
		duplicate/d/o lambda, $(name+"_L")
		duplicate/d/o timeData, $(name+"_T")
		duplicate/d/o spec $(name+"_raw")
		WAVE spec=$(name+"_raw"), spec2=$(name), lambda=$(name+"_L"), timeData=$(name+"_T")
		Killwaves $("CP_"+name), $("RP_"+name)
		WAVE spec=$(name), lambda=$("CP_"+name), timeData=$("RP_"+name)
		WAVE spec=$(name+"_raw"), spec2=$(name), lambda=$(name+"_L"), timeData=$(name+"_T")
	
		//Make eV Wave/
		Duplicate/D/O lambda,eV
		Wave eV
		eV = 1250/lambda
		DeletePoints/M=0 0,6,eV
		DeletePoints/M=0 774,1000,eV
		Duplicate/D/O eV,shortWL
		Wave shortWL
		shortWL = 1250/eV

		//Scale data for time
		Duplicate/d/o timeData dt
		Differentiate dt
		wavestats/q dt
		print "Average Data Collection Rate(ms)=", V_avg*1000 //convert to seconds from ms
		variable avgDt=V_avg 
		setscale/p x, 0, avgDt, spec, spec2
	
		//interpolate data for wavelength
		Duplicate/d/o lambda dummyspecRaw, dummyspecScaled
		variable firstLambda=lambda[0]
		variable lastLambda=lambda[numpnts(lambda)-1]
		setscale/I y, firstLambda, lastLambda, spec2
		setscale/I x, firstLambda, lastLambda, dummyspecScaled
		Variable i
		For(i=0; i<numpnts(timedata); i+=1 )
			multithread dummyspecRaw=spec[i][p]
			multithread dummyspecScaled=interp(x, lambda, dummyspecRaw)
			multithread spec2[i][]=dummyspecScaled[q]
			If( i>0 && mod(i,2000)==0 )
				print "Processing time index="+num2str(i)
			endif
		endfor
		Killwaves dummyspecRaw, dummyspecScaled, dt,spec
	
		Wave wave2D = $name
		Variable tf = DimSize(wave2D,0)	//What is the last spectra collected?
		If( d )
			dowindow/F $name
			If(!V_Flag)
				Display /W=(10,50,400,340) /N=$name as name
				appendimage spec2
				ModifyGraph grid=2,tick=2,minor=1,standoff=0
				Label left "Wavelength [nm]"
				Label bottom "Time [sec]"
				ModifyGraph margin(left)=40,margin(bottom)=33,margin(top)=14,margin(right)=70
				ColorScale/C/N=text0 fsize=10,fstyle=1, "Absorbance[a.u.]"	
				ColorScale/C/N=text0/Z=0/A=RC/X=-24/Y=-5.02	
				WaveStats/Q/RMD=[tf,tf][65,773] wave2D	//Finds max absorbance b/n 350-900nm for last collected spectra
				Variable maxAbs = V_max
				ModifyImage $name ctab={0,maxAbs,Terrain,0}
				SetAxis left wli,wlf	
			endif
		endif
		
		//Extract the individual absorbance profiles from 2D wave
		Variable tms = V_avg*1000 //Data collection rate in ms
		Variable ti
		Prompt ti, "Enter starting point for data analysis: "
		DoPrompt "Enter starting point",ti
		if (V_Flag)
			return -1								// User canceled
		endif
		Wave wave2D = $name
		tf = DimSize(wave2D,0)
		
		//Check that no processed waves have been created yet
		String oldSpecs = WaveList(name+"_Prof*",";","")
		Variable nOldSpec = ItemsInList(oldSpecs)
		for(i=0;i<nOldSpec;i+=1)
			String cSpecName = StringFromList(i,oldSpecs)
			KillWaves $cSpecName
		endfor
		ExtractAll(wave2D,tms,ti,tf)

		String specList = WaveList(name+"_Prof*",";","")
		Variable/G nSpec =ItemsInList(specList)-1

		Wave twave,SigCs		
		processTRabsorbance(name,eV,twave,nSpec,name)
		print "Your data has been processed"
	endif
	
	if(doFits)
		SetDataFolder $fname
		String path2nSpec = fname +":nSpec"
		NVAR nSpec = $path2nSpec 
		//What kind of aggregate is present? H or J?
		if(StringMatch(aggregateType,"H"))
			Wave twave,SigCs			
			SpanoAnalysis(name,eV,tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale,twave,SigCs,nSpec)
			print "Data has been fitted to Spano Model"
			print "Your data has been processed"
		elseif(StringMatch(aggregateType,"J"))
			print "This feature has not been implemented yet"
		endif
	endif
							
	SetDataFolder $CurrentFolder
End

Function ExtractAll(image,tf,ini,final)
	
	Wave image
	Variable tf,ini,final
	Variable S,W,Ep,E00,sig,m,scale
	
	//Generates list of values with between ini and final in increments of 1
	Variable i
	
	String test, longtest = ""
	
	Variable dif = final - ini
	//How many spectra to process? If there are more than 1000 spectra subsequent analysis
	//can crash IGOR...
	if(dif < 500)
		for(i=(ini);i<=(final);i=i+1) 
			test= num2str(i)
			longtest = addlistitem(test, longtest, ";", inf)
		endfor	
	elseif(dif > 500)
		for(i=(ini);i<=(final);i+=16) 
			test= num2str(i)
			longtest = addlistitem(test, longtest, ";", inf)
		endfor	
	endif 	
	
	Wave spec0
	
	Make/O/N=1044 spec0					//Wave with 1044 rows b/c of wavelength
	
	for(i=0;i<ItemsInList(longtest);i=i+1)
		spec0 = image[str2num(stringfromlist(i,longtest,";"))][p]
			String  newspec=NameOfWave(image)+"_Prof"+num2str(i)
			Duplicate/O spec0, $(newspec)
	endfor
	KillWaves spec0
	
	//Make timewave for analysis of film dynamics
	Make/O/N=(ItemsInList(longtest)) rawtime
	Wave rawtime
	for(i=0;i<ItemsInList(longtest);i=i+1)
		rawtime[i]=str2num(stringfromlist(i,longtest,";"))
	endfor
	
	Make/O/N=(ItemsInList(longtest)) twave
	Wave twave
	twave= (rawtime-rawtime[0])*tf

	//Make coefficient wave for Sigmoid Fits
	Make/D/O/N=4 SigCs
	SigCs = {10,10,10,10}
	
	KillWaves rawtime
End


Function SpanoAnalysis(w0,xwave,tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale,twave,sigC,fnum) //Carries out the entire Spano Analysis

	String w0							//w0 is the basename of the absorbance profiles
	Wave xwave,twave,sigC// xwave is the trimmed eV wave, cwave is the Spano Cs, twave is the time wave, Sigcwave is the Sigmoid Cs
	Variable tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale,fnum						//last number of absorbance profile to analyze
	
	//This loop subtracts background, normalizes, fits Spano model, and determines amorphous spectra for all collected spectra
	//It also plots each of the aforementioned data into respective plots
	Variable i
	String rawAbs, normAbs,spanoAbs,amorphousAbs
	//Process the data
	for(i=0; i<=fnum; i=i+1) 			
		normAbs = w0 + "_Prof" + num2str(i) + "_n"	//Normalized absorbance spectra
		Wave w2 = $normAbs 
		SpanoFit(w2,xwave,tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale)
		spanoAbs = "Spano_" + w0 + "_Prof" + num2str(i) + "_n"
		Wave w3 = $spanoAbs
		AmorphousSpectra(w2,w3)
		amorphousAbs =  w0 + "_Prof" + num2str(i)+"_n_A"
		Wave w4 = $amorphousAbs	
	endfor
	//Plot the spectra
	Variable d=1
	if(d)
		String spanoPltName = w0 + "_SpanoFits"
		DoWindow $spanoPltName
		if(!V_Flag)
		 	Display/N=$spanoPltName/W=(0,0,600,300)
			for(i=0;i<=fnum;i+=1)
				spanoAbs = "Spano_" + w0 + "_Prof" + num2str(i) + "_n"
				Wave w3 =$spanoAbs
				AppendToGraph/W=$spanoPltName w3 vs xwave
				ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1
				Label bottom "eV"
				Label left "Absorbance"
				SetAxis left 0,1.2
				ApplyColorTableToTopGraphSpano("Rainbow")
				ColorScale/C/N=text0/A=LT  ctab={0,100,Rainbow,0},lblRot=180
				ColorScale/C/N=text0 "Drying Time -->"
			endfor
		endif
		
		String amorphPltName = w0 + "AmorphousSpec"
		DoWindow $amorphPltName
		if(!V_Flag)
		 	Display/N=$amorphPltName/W=(0,0,600,300)
			for(i=0;i<=fnum;i+=1)
				amorphousAbs =  w0 + "_Prof" + num2str(i) +"_n_A"
				Wave w4 =$amorphousAbs
				AppendToGraph/W=$amorphPltName w4 vs xwave
				ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1
				Label bottom "eV"
				Label left "Absorbance"
				SetAxis left 0,1.2
				ApplyColorTableToTopGraphSpano("Rainbow")
				ColorScale/C/N=text0/A=LT  ctab={0,100,Rainbow,0},lblRot=180
				ColorScale/C/N=text0 "Drying Time -->"
			endfor
		endif
	endif
	
	//Determine Oscillator Strengths
	String osTotal = "osTotal_" + w0
	String osAggregate = "osAgg_" + w0 
	String osAmorphous = "osAmo_" + w0
	String pAgg = "pAgg_" + w0
	String pAmo = "pAmo_" + w0
	Make/O/N=(fnum+1) $osTotal, $osAggregate, $osAmorphous,$pAgg, $pAmo
	Wave osTot = $osTotal,osAgg = $osAggregate,osAmo =$osAmorphous,perAgg=$pAgg,perAmo=$pAmo
	
	for(i=0;i<=fnum;i+=1)
		normAbs = w0 + "_Prof" + num2str(i) + "_n"
		spanoAbs = "Spano_" + w0 + "_Prof" + num2str(i) + "_n"
		amorphousAbs =  w0 + "_Prof" + num2str(i) +"_n_A"
		Wave w2 =$normAbs,w3 =$spanoAbs, w4 =$amorphousAbs
		osTot[i]=Area(w2)	
		osAgg[i]=Area(w3) 
		osAmo[i]=Area(w4)	
	endfor
	
	//Calculate and plot percent aggregate/amorphous chains in film
	perAgg = (osAgg/osTot)*100
	perAmo = (osAmo/osTot)*100
	
	//Extract and Graph Spano parameters
	String HRF = "S_" + w0
	String EB = "W_" + w0
	String GW = "ED_" + w0
	String ZeroVibE = "E00_"+ w0
	String Amplitude = "Amp_" + w0
	Make/O/N=(fnum+1) $HRF,$EB,$GW,$ZeroVibE,$Amplitude
	Wave S = $HRF,EBW = $EB, ED =$GW, E_00 = $ZeroVibE, Amp = $Amplitude
	
	for(i=0;i<=fnum;i+=1)
		String SpanoCwName = "Cs_"+ w0 + "_Prof"+ num2str(i) + "_n"
		Wave spanoCwave = $SpanoCwName
		S[i]   = spanoCwave[0]
		EBW[i]   = spanoCwave[1]
		E_00[i] = spanoCwave[3]
		ED[i]  = spanoCwave[4]
		Amp[i] = spanoCwave[6]	
	endfor
	
	Variable difPnts = abs(numPnts(twave) - numpnts(S))	
	DeletePoints	numpnts(twave),difPnts,S,EBW,ED,E_00,Amp
	
	if(d)
		//Plot results of Spano Fits
		String sPltName = w0 + "_HRF"
 		Display/N=$sPltName/W=(0,0,200,200) S vs twave 
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,mode=3,marker=19
		Label bottom "Time[ms]\U"
		Label left "S"
		
		String wPltName = w0 + "_ExcitonBandWidth"
	 	Display/N=$wPltName/W=(200,0,400,200) EBW vs twave
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,mode=3,marker=19
		Label bottom "Time[ms]\U"
		Label left "W[meV]"
		
		String e00PltName = w0 + "_E00"
 		Display/N=$e00PltName/W=(400,0,600,200) E_00 vs twave
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,mode=3,marker=19
		Label bottom "Time[ms]\U"
		Label left "E\\B0-0\\M[eV]"
		
		String disPltName = w0 + "_EnergeticDisorder"
		Display/N=$disPltName/W=(600,0,800,200) ED vs twave
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,mode=3,marker=19
		Label bottom "Time[ms]\U"
		Label left "σ[meV]"
		
		String ampPltName = w0 + "_ProportionalityFactor"
	 	Display/N=$ampPltName/W=(800,0,1000,200) Amp vs twave
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,mode=3,marker=19
		Label bottom "Time[ms]\U"
		Label left "Amplitude"
		
		String pAggPltName = w0 + "_PercentAggregates"
		Display/N=$pAggPltName/W=(0,200,200,400)
		AppendToGraph perAgg,perAmo vs twave
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,mode=3,marker=8, rgb($pAmo)=(0,0,0),rgb($pAgg)=(39321,1,1)
		Label bottom "Time[ms]\U"
		Label left "%"	
	endif
End

Function processTRabsorbance(w0,xwave,twave,fnum,name) //Processes time resolved absorbance data
	String w0							//w0 is the basename of the absorbance profiles
	Wave xwave,twave// xwave is the trimmed eV wave,twave is the time wave
	Variable fnum						//last number of absorbance profile to analyze
	String name
	//This loop subtracts background, normalizes, fits Spano model, and determines amorphous spectra for all collected spectra
	//It also plots each of the aforementioned data into respective plots
	Variable i
	String rawAbs, normAbs
	//Process the data
	for(i=0; i<=fnum; i=i+1) 			
		rawAbs = w0 + "_Prof" + num2str(i)
		Wave w1=$rawAbs
		BaselineSubtraction(w1)
		normAbs = w0 + "_Prof" + num2str(i) + "_n"	//Normalized absorbance spectra
		Wave w2 = $normAbs 	
	endfor
	//Plot the spectra
	Variable d=1
	if(d)
	String rawSpecPltName = name + "_RawSpec"
	DoWindow/F $rawSpecPltName
	if(V_Flag!=1)
	 	Display/N=$rawSpecPltName/W=(0,0,600,300) 
		for(i=0;i<=fnum;i+=1)
			rawAbs = w0 + "_Prof" + num2str(i)
			Wave w1 =$rawAbs
			AppendToGraph/W=$rawSpecPltName w1 vs xwave
			ModifyGraph/W=$rawSpecPltName mirror=1,nticks=10,minor=1,fStyle=1
			Label/W=$rawSpecPltName bottom "Photon Energy[eV]"
			Label/W=$rawSpecPltName left "Absorbancee[a.u.]"
			SetAxis/W=$rawSpecPltName left 0,1.2
			ApplyColorTableToTopGraphSpano("Rainbow")
			ColorScale/C/N=text0/A=LT  ctab={0,100,Rainbow,0},lblRot=180
			ColorScale/C/N=text0 "Drying Time -->"	
		endfor
	endif
	
	String normSpecPltName = name + "_NormSpec"
	DoWindow/F $normSpecPltName
	if(V_Flag!=1)
	 	Display/N=$normSpecPltName/W=(600,0,1200,300) 
		for(i=0;i<=fnum;i+=1)
			normAbs = w0 + "_Prof" + num2str(i) + "_n"
			Wave w2 =$normAbs
			AppendToGraph/W=$normSpecPltName w2 vs xwave
			ModifyGraph/W=$normSpecPltName mirror=1,nticks=10,minor=1,fStyle=1
			Label/W=$normSpecPltName bottom "Photon Energy[eV]"
			Label/W=$normSpecPltName left "Absorbance[a.u.]"
			SetAxis/W=$normSpecPltName left 0,1.2
			ApplyColorTableToTopGraphSpano("Rainbow")
			ColorScale/C/N=text0/A=LT  ctab={0,100,Rainbow,0},lblRot=180
			ColorScale/C/N=text0 "Drying Time -->"	
		endfor
	endif
	endif
End

Function ApplyColorTableToTopGraphSpano(ctabname)
    String ctabname

    String graphName = WinName(0, 1)
    if (strlen(graphName) == 0)
        return -1
    endif

    Variable numTraces = ItemsInList(TraceNameList(graphName,";",3))

    if (numTraces <= 0)
        return -1
    endif
   
    Variable denominator= numTraces-1
    if( denominator < 1 )
        denominator= 1    // avoid divide by zero, use just the first color for 1 trace
    endif

    ColorTab2Wave $ctabname // creates M_colors
    Wave M_colors
    Variable numRows= DimSize(M_colors,0)
    Variable red, green, blue
    Variable i, index
    for(i=0; i<numTraces; i+=1)
        index = round(i/denominator * (numRows-1))  // spread entire color range over all traces.
        ModifyGraph/W=$graphName rgb[i]=(M_colors[index][0], M_colors[index][1], M_colors[index][2])
    endfor
    return 0
End


//Functions referenced during SpanoAnalysis

Function BaselineSubtraction(w) 	//This function finds the minimum value over a specified range and subtracts it from the parent wave. 
															// The baseline corrected wave is then smoothed out.
	WAVE w												// w1 is the initial wave to be baselined
	
	Variable leftmark=65				// leftmark is the leftmost "x" value of the range, rightmark is the rightmost "x" value of the range
	Variable rightmark=800			//rangetype should be a point in the wave that is within the specified range
	Variable rangetype=65
	
	String blName= NameOfWave(w) +"_bl" 			//Renames the original wave 
	String nName = NameOfWave(w) +"_n"
	
	if (rangetype<65 || rangetype>800)				//if the value of for rangetype is outside the specified range, return null result
		return -1									//The range is setup b/n 300 and 629 because that is the corresponding ev/nm range that I care about in the 
	endif										// new waves
	
	Duplicate/O w $blName						//Duplicates original wave and renames it according to the outputName function
	Wave bl = $blName
	
	if (rangetype==0)								//Determines the various wave statistics over the specified range
		WaveStats/Q bl
	endif
	
	if (rangetype==65)
		WaveStats/Q/R=[leftmark,rightmark] bl
	endif
	
	if (rangetype==800)
		WaveStats/Q/R=(leftmark,rightmark) bl
	endif
	
	bl = bl - V_min							//Carry out baseline subtraction by utlizing V_min determined from previous loops
	DeletePoints/M=0 0,6,bl					//Remove irrelevant/uninteresting points from polished data
	DeletePoints/M=0 774,1000,bl
	//Smooth 10,output						//Smoothes out the baseline corrected wave. 
	
	Duplicate/O bl $nName
	Wave n = $nName
	
	WaveStats/Q n
	n = bl/V_max
End

Function AmorphousSpectra(w1,w2)		//Determines amorphous spectra by subtracting Spano fit from total spectra

	WAVE w1,w2						//w1=TotalSpectra and w2=Spano Fit
	String outputName= NameOfWave(w1)+"_A"
	
	Duplicate/O w1 $(outputName)
	Wave output = $(outputName)

	output = w1 - w2	
End

Function SpanoFit(w1,xwave,tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale)

	Wave w1,xwave
	Variable tol,HRS,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale
	
	String outputName="Spano_" + NameOfWave(w1) 
	String cwaveName="Cs_" +NameOfWave(w1)
	String sigmaName = "Sigma_"+NameOfWave(w1)
	String resName = "res_"+NameOfWave(w1)
	//Variable/G V_FitTol=0.0001
	//Variable/G V_FitMxIters=150
	//Prepare parameter wave,constraints, hold string, and epsilon wave
	Wave T_ConstraintsSpano = makeSpanoConstraints(S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale)
	Wave cwave = makeSpanoPWave(HRS,W,Ep,E00,sig,m,scale)
	Wave EpsilonSpano = makeSpanoEpsWave(tol,holdS,holdW,holdE00,holdSig,holdScale)
	Duplicate /O w1, $(outputName), $(resName)
	Wave output=$(outputName)
	Wave res = $(resName)
	output = 0
	String H = makeSpanoHoldString(holdS,holdW,holdE00,holdSig,holdScale)
	//Do the fit
	Variable/G V_FitError = 0
	FuncFit/L=624/Q /H=H  Spano cwave w1[317,435] /X=xwave/E=EpsilonSpano/C=T_ConstraintsSpano/D=output/R=res
	output= Spano(cwave,xwave)
	
	Duplicate/O cwave $(cwaveName)
	Wave output1=$(cwaveName)
	
//	Wave W_Sigma
//	Duplicate/O W_Sigma $(sigmaName)
//	Wave Uncertainties = $(sigmaName)

	//SpanoGaussians("P3HT",output1,xWave)
End

Function SpanoFit2(w1,xwave,cWave,tol)

	Wave w1,xwave,cwave
	Variable tol
	
	String outputName="Spano_" + NameOfWave(w1) 
	String cwaveName="Cs_" +NameOfWave(w1)
	String sigmaName = "Sigma_"+NameOfWave(w1)
	
	Variable/G V_FitTol=0.000001
	Variable/G V_FitMxIters=150
	
	Make/D/O/N=7 EpsilonSpano
	EpsilonSpano[0]={0,10^-tol,0,10^-tol,10^-tol,10^-tol,10^-tol}
	
	Make/D/O/T/N=7 T_ConstraintsSpano
	T_ConstraintsSpano[0]={"K1 > 0.04","K1 < 0.1","K3 > 1.85","K3 < 2.2","K4 > 0.07","K4 < 0.09","K6 > 0","K6 < 3"}
	Duplicate /O w1 $(outputName)
	Wave output=$(outputName)
	
	FuncFit/Q/N=1 /L=624 /H="1010010"  Spano cwave w1[317,435] /X=xwave /D/E=EpsilonSpano/C=T_ConstraintsSpano//F={0.95,7}/R 	//[600,775] for S-Flame [1000,1225] for T-Flame
	output= Spano(cwave,xwave)
	
	Duplicate/O cwave $(cwaveName)
	Wave output1=$(cwaveName)
	
	Wave W_Sigma
	Duplicate/O W_Sigma $(sigmaName)
	Wave Uncertainties = $(sigmaName)

	//SpanoGaussians("P3HT",output1,xWave)
End

Function Spano(w,E) : FitFunc
	Wave w
	Variable E
	Variable m,n
	Variable Total
	Variable Amp,Gterm,SumTerm,Form
	//w[0] == Huang Rhys Factor //w[1] == Exciton Bandwidth
	//w[2] == Intermolecular Vibrational Energy //w[3] ==Zero-zero Vibrational Transition Energy
	//w[4] == Energetic Disorder //w[5] == Number of Vibrational Levels //w[6] == Proportionality Factor
	for(m=0;m<=w[5];m+=1)
			//Transition can't occur in same vibrational level										
				Amp = w[6]*(exp(-w[0])*(w[0]^m)/(factorial(m)))
				GTerm = (w[1]*Exp(-w[0])/(2*w[2]))
				n=0
				SumTerm = 0
				//This term needs to be within its own loop to prevent overaddition
				for(n=0;n<=w[5];n+=1)	
					if (m != n)
						SumTerm += (w[0]^n)/(factorial(n)*(n-m))
					else
						SumTerm += 0
					endif
				endfor
				Form = exp(-(((E-w[3]-m*w[2]-(1/2)*w[1]*((w[0])^m)*exp(-w[0]))^2)/(2*(w[4])^2)))
				
				total += Amp* (( 1-Gterm * SumTerm)^2)*Form
	endfor
	return total		
End

Function SpanoGaussians(Name,w,E)
	String Name
	Wave w,E
	
	variable m,n
	Make/N=774/O Amp,Gterm,SumTerm,Form
	
	String CurrentFOlder = GetDataFOlder(1)
	//NewDataFolder/o/s $"root:Spano_Gaussians"
	//NewDataFolder/o/s $name
	
	Make/N=774/o Total_Gaussian=0

	String GaussName
	for(m=0;m<=w[5];m+=1)
		GaussName = Name + "_Gauss"+"_"+num2str(m)
		Make/o/N=774 $GaussName
	endfor
	
	//Make Gaussian Components
		for(m=0;m<=w[5];m+=1)
			//Transition can't occur in same vibrational level						
				String CurrentGauss = Name + "_Gauss"+"_"+num2str(m)
				wave Gaussian = $CurrentGauss
				
				Amp = w[6]*(exp(-w[0])*(w[0]^m)/(factorial(m)))
				GTerm = (w[1]*Exp(-w[0])/(2*w[2]))
				n=0
				SumTerm = 0
				for(n=0;n<=w[5];n+=1)
					if (m != n)
						SumTerm += (w[0]^n)/(factorial(n)*(n-m))
					else
						SumTerm += 0
					endif
				endfor
				Form = exp(-(((E-w[3]-m*w[2]-(1/2)*w[1]*((w[0])^m)*exp(-w[0]))^2)/(2*(w[4])^2)))
				
				total_Gaussian += Amp* (( 1-Gterm * SumTerm)^2)*Form
				Gaussian =  Amp* (( 1-Gterm * SumTerm)^2)*Form
			//	Appendtograph Gaussian vs E
			//	ModifyGraph rgb($Nameofwave(Gaussian))=(0,0,0)
	endfor
	//Appendtograph Total_Gaussian vs E
	
	//SetDataFolder $CurrentFolder
end

Function GetW(yw,xw,aggType)

	Wave yw,xw	//yw is the absorbance spectra to process. xw is the wavelength in nm
	String aggType
	
	Variable W_H,W_J
	Variable Ep = 0.179, A00wl=500,A01wl=600
	
	Prompt A00wl, "What is the wavelength of the 0-0 transition?"
	Prompt A01wl, "What is the wavelength of the 0-1 transition?"
	DoPrompt "Which aggregate type is present? Enter H or J:",A00wl,A01wl
	//Find the corresponding points for the wavelength of the A00 and A01 transition
	Variable pA00 = BinarySearch(xw,A00wl)
	Variable pA01 = BinarySearch(xw,A01wl)
	//Use the point values from above to determine the corresponding absorbance intensity
	Variable A00 = yw[pA00]
	Variable A01 = yw[pA01]
	if(stringMatch(aggtype, "H"))
		W_H = (Ep*(sqrt((A00/A01))-1))/(0.073*sqrt((A00/A01))-0.24)
		PRINT W_H 
	elseif(stringMatch(aggtype, "J"))
		W_J = -((Ep*(sqrt((A00/A01))-1))/(0.073*sqrt((A00/A01))+0.24))
		PRINT W_J
	else
		return -1
	endif
	
	Wave W
	if(!WaveExists(W))
		Make/O/N=1 W
	endif		
End

Menu "Macros"
	"Spano Panel", SpanoPanel()
end

Function SpanoPanel()

	String panelName = "Spano_Panel"
	DoWindow $panelName
	if(V_Flag)
		KillWindow $panelName
	endif
	NewPanel/K=1/N=$panelName/W=(0,0,400,350)
	
	String fileName = ""
	Variable HRS   = 1
	Variable W   = 0.05
	Variable Ep  = 0.179
	Variable E00 = 2.00
	Variable sig = 0.05
	Variable m   = 5
	Variable amp = 1
	Variable tol = 1E-7
	Variable S_lo = 0.5
	Variable S_hi = 3
	Variable W_lo = 0.0001
	Variable W_hi = 2
	Variable E00_lo = 1.7
	Variable E00_hi = 3.8
	Variable sig_lo = 0.0001
	Variable sig_hi = 2
	Variable amp_lo = 0.001
	Variable amp_hi = 5
	//Data loading section
	TitleBox DataLoading        title="Data Loading"                    ,pos={169,4}
	SetVariable fileName        title="File Name"         ,size={175,20},pos={5,33} ,limits={-inf,inf,0},value=_STR:fileName
	CheckBox display2D title="Display 2D Absorbance?"                   ,pos={198,33},value=1
	
	//Fit parameter section
	TitleBox SpanoParams title="Spano Parameters"  ,pos={159,57}
	SetVariable HRFactor         title="S"         ,pos={59,105} ,size={50,20} ,limits={0,3,0.1}      ,value=_NUM:HRS
	SetVariable excitonBandwidth title="W[eV]"     ,pos={29,127} ,size={80,20} ,limits={0,1,0.01}     ,value=_NUM:W
	SetVariable IMVibE           title="Ep[eV]"    ,pos={29,148} ,size={80,20} ,limits={-inf,inf,0}   ,value=_NUM:Ep
	SetVariable E0_0             title="E\B0-0[eV]"  ,pos={10,169} ,size={100,20},limits={1.8,3.1,0.01} ,value=_NUM:E00
	SetVariable disorder         title="σ[eV]"     ,pos={30,191} ,size={80,20} ,limits={0,1,0.01}     ,value=_NUM:sig
	SetVariable nGauss           title="m"         ,pos={60,212} ,size={50,20} ,limits={1,10,1}       ,value=_NUM:m
	SetVariable scale            title="Scaling"   ,pos={11,233} ,size={100,20},limits={0,10,1}       ,value=_NUM:amp	
	SetVariable tol              title="Tolerance" ,pos={160,82} ,size={150,20}  ,limits={1E-15,0.01,1E-12}       ,value=_NUM:tol
	
	CheckBox hS title="Hold S?"         ,pos={121,105},value=1
	CheckBox hW title="Hold W?"         ,pos={121,127},value=0
	CheckBox hE title="Hold E\B0-0?"    ,pos={121,169},value=0
	CheckBox hD title="Hold σ?"         ,pos={121,191},value=0
	CheckBox hA title="Hold Scale?"     ,pos={121,233},value=0
	
	SetVariable S_lo         title="S\BLO"      ,pos={213,105}  ,size={75,20}  ,limits={0,3,0.1}          ,value=_NUM:S_lo
	SetVariable S_hi         title="S\BHI"      ,pos={296,105}  ,size={75,20}  ,limits={0,1,0.01}         ,value=_NUM:S_hi
	SetVariable W_lo         title="W\BLO"      ,pos={213,127}  ,size={75,20} ,limits={-inf,inf,0}       ,value=_NUM:W_lo
	SetVariable W_hi         title="W\BHI"      ,pos={296,127}  ,size={75,20} ,limits={1.8,3.1,0.01}      ,value=_NUM:W_hi
	SetVariable E00_lo       title="E\B0-0,LO"  ,pos={213,169}  ,size={75,20} ,limits={0,1,0.01}        ,value=_NUM:E00_lo
	SetVariable E00_hi       title="E\B0-0,HI"  ,pos={296,169}  ,size={75,20} ,limits={1,10,1}          ,value=_NUM:E00_hi
	SetVariable sig_lo       title="σ\BLO"      ,pos={213,191}  ,size={75,20} ,limits={0,10,1}            ,value=_NUM:sig_lo	
	SetVariable sig_hi       title="σ\BHI"      ,pos={296,191} ,size={75,20} ,limits={1E-15,0.01,1E-12}  ,value=_NUM:sig_hi
	SetVariable A_lo         title="A\BLO"      ,pos={213,233}  ,size={75,20} ,limits={0,3,0.1}         ,value=_NUM:amp_lo
	SetVariable A_hi         title="A\BHI"      ,pos={296,233}  ,size={75,20} ,limits={0,1,0.01}        ,value=_NUM:amp_hi
	
	CheckBox analyze title="Run Spano Analysis?"         ,pos={144,270},value=1
	CheckBox load    title="Load 2D Absorbance?"         ,pos={5,270},value=1
	Button loadTimeSeriesButton title="Run",size={125,20},pos={272,266},proc=load2DAbsBut
	
	//Plot Spano Results
	TitleBox pltTitle   title="Plot Spano Set"  ,pos={165,289}
	SetVariable setVal  title="Spec #"        ,pos={82,321}  ,size={100,20} ,limits={0,inf,1}        ,value=_NUM:0,proc=SetVarSpano1
	Button pltSet       title="Plot Spano Set",size={125,20},pos={195,319},proc=pltSpanoSetBut
End

Function load2DAbsBut(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			String panelName = "Spano_Panel"
			ControlInfo/W=$panelName fileName
			String name = S_Value		
			ControlInfo/W=$panelName tol
			Variable tol = V_Value			
			ControlInfo/W=$panelName display2D
			Variable display2D = V_Value			
			ControlInfo/W=$panelName HRFactor
			Variable HRFactor = V_Value			
			ControlInfo/W=$panelName excitonBandwidth
			Variable W = V_Value			
			ControlInfo/W=$panelName IMVibE
			Variable Ep = V_Value			
			ControlInfo/W=$panelName E0_0
			Variable E00 = V_Value			
			ControlInfo/W=$panelName disorder
			Variable sig = V_Value			
			ControlInfo/W=$panelName nGauss
			Variable m = V_Value			
			ControlInfo/W=$panelName scale
			Variable scale = V_Value
			ControlInfo/W=$panelName S_lo
			Variable S_lo = V_Value
			ControlInfo/W=$panelName S_hi
			Variable S_hi = V_Value
			ControlInfo/W=$panelName W_lo
			Variable W_lo = V_Value
			ControlInfo/W=$panelName W_hi
			Variable W_hi = V_Value
			ControlInfo/W=$panelName E00_lo
			Variable E00_lo = V_Value
			ControlInfo/W=$panelName E00_hi
			Variable E00_hi = V_Value
			ControlInfo/W=$panelName sig_lo
			Variable sig_lo = V_Value
			ControlInfo/W=$panelName sig_hi
			Variable sig_hi = V_Value
			ControlInfo/W=$panelName A_lo
			Variable A_lo = V_Value
			ControlInfo/W=$panelName A_hi
			Variable A_hi = V_Value
			ControlInfo/W=$panelName hS
			Variable holdS = V_Value
			ControlInfo/W=$panelName hW
			Variable holdW = V_Value
			ControlInfo/W=$panelName hE
			Variable holdE00 = V_Value
			ControlInfo/W=$panelName hD
			Variable holdSig = V_Value
			ControlInfo/W=$panelName hA
			Variable holdScale = V_Value
			ControlInfo/W=$panelName analyze
			Variable analyze = V_Value
			ControlInfo/W=$panelName load
			Variable load = V_Value
			ITS(name,tol,HRFactor,W,Ep,E00,sig,m,scale,S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale,d=display2D,load=load,doFits=analyze)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function/WAVE makeSpanoConstraints(S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale)

	Variable S_lo,S_hi,W_lo,W_hi,E00_lo,E00_hi,sig_lo,sig_hi,A_lo,A_hi,holdS,holdW,holdE00,holdSig,holdScale
	
	Make/O/T/N=10 SpanoConstraints
	SpanoConstraints[0] = "K0 >= " + num2str(S_lo)
	SpanoConstraints[1] = "K0 <= " + num2str(S_hi)
	SpanoConstraints[2] = "K1 >= " + num2str(W_lo)
	SpanoConstraints[3] = "K1 <= " + num2str(W_hi)
	SpanoConstraints[4] = "K3 >= " + num2str(E00_lo)
	SpanoConstraints[5] = "K3 <= " + num2str(E00_hi)
	SpanoConstraints[6] = "K4 >= " + num2str(sig_lo)
	SpanoConstraints[7] = "K4 <= " + num2str(sig_hi)
	SpanoConstraints[8] = "K6 >= " + num2str(A_lo)
	SpanoConstraints[9] = "K6 <= " + num2str(A_hi)
	
	if(holdScale)
		DeletePoints 8,2,SpanoConstraints
	endif
	
	if(holdSig)
		DeletePoints 6,2,SpanoConstraints
	endif
	
	if(holdE00)
		DeletePoints 4,2,SpanoConstraints
	endif
	
	if(holdW)
		DeletePoints 2,2,SpanoConstraints
	endif
	
	if(holdS)
		DeletePoints 0,2,SpanoConstraints
	endif
	
	return SpanoConstraints

End

Function/WAVE makeSpanoPWave(HRS,W,Ep,E00,sig,m,scale)

	Variable HRS,W,Ep,E00,sig,m,scale
	
	Make/O/D/N=7 SpanoWave
	SpanoWave[0] = HRS
	SpanoWave[1] = W
	SpanoWave[2] = Ep
	SpanoWave[3] = E00
	SpanoWave[4] = sig
	SpanoWave[5] = m
	SpanoWave[6] = scale
	
	return SpanoWave
End

Function/S makeSpanoHoldString(holdS,holdW,holdE00,holdSig,holdScale)

	Variable holdS,holdW,holdE00,holdSig,holdScale
	
	String H = ""
	
	//Hold Huang Rhys Factor?
	if(holdS)
		H+= "1"
	else
		H+= "0"
	endif
	
	//Hold Exciton Bandwidth?
	if(holdW)
		H+= "1"
	else
		H+= "0"
	endif
	
	//Intermolecular vibrational energy is always held
	H+="1"
	
	//Hold Zeroth Vibrational Energy?
	if(holdE00)
		H+= "1"
	else
		H+= "0"
	endif
	
	//Hold energetic disorder?
	if(holdSig)
		H+= "1"
	else
		H+= "0"
	endif
	
	//Number of Gaussians to fit is always held
	H+="1"
	
	//Hold scaling factor?
	if(holdScale)
		H+= "1"
	else
		H+= "0"
	endif
	
	return H

End

Function/WAVE makeSpanoEpsWave(tol,holdS,holdW,holdE00,holdSig,holdScale)

	Variable tol,holdS,holdW,holdE00,holdSig,holdScale
	
	Make/O/N=7 eps
	
	//Hold Huang Rhys Factor?
	if(holdS)
		eps[0] = 0
	else
		eps[0] = tol
	endif
	
	//Hold Exciton Bandwidth?
	if(holdW)
		eps[1] = 0
	else
		eps[1] = tol
	endif
	
	//Intermolecular vibrational energy is always held
	eps[2] = 0
	
	//Hold Zeroth Vibrational Energy?
	if(holdE00)
		eps[3] = 0
	else
		eps[3] = tol
	endif
	
	//Hold energetic disorder?
	if(holdSig)
		eps[4] = 0
	else
		eps[4] = tol
	endif
	
	//Number of Gaussians to fit is always held
	eps[5] = 0
	
	//Hold scaling factor?
	if(holdScale)
		eps[6] = 0
	else
		eps[6] = tol
	endif
	
	return eps

End

Function plotSpanoSet(name,val)
	
	String name
	Variable val
	
	String fName = "root:" + name +":"
	SetDataFolder root:
	SetDataFolder $fName
	String amoName = name + "_Prof" + num2str(val) + "_n_A"
	String spaName = "Spano_" + name + "_Prof" + num2str(val) + "_n"
	String totName = name + "_Prof" + num2str(val) + "_n"
	String cofName = "Cs_" + name + "_Prof" + num2str(val) + "_n"
	Wave amo = $amoName
	Wave spa = $spaName
	Wave tot = $totName
	Wave cof = $cofName
	Wave eV
	Variable i
	//Make Gaussians for current set
	SpanoGaussians(name,cof,eV)
	String gaussList = WaveList( Name + "_Gauss*",";","")
	Variable nGauss = ItemsInList(gaussList)
	String pltName = "Spano_Summary" 
	DoWindow	$pltName
	if(!V_Flag)
		Display/N=$pltName/K=1 amo,spa,tot vs eV
		ModifyGraph tick=2,mirror=1,minor=1,fStyle=1
		Label left "Absorbance [a.u.]"
		Label bottom "Photon Energy [eV]"
		SetAxis bottom 1.6,3.5
		ModifyGraph lsize=1.5,rgb($spaName)=(1,26221,39321),rgb($totName)=(0,0,0)
		for(i=0;i<nGauss;i+=1)
			String cGauss = StringFromList(i,gaussList)
			Wave g = $cGauss
			AppendToGraph g vs eV
			ModifyGraph lstyle($cGauss)=3,rgb($cGauss)=(21845,21845,21845),lsize($cGauss)=1.5
		endfor
		Legend/C/N=text0/J/A=RT "\\JCSet # "+num2str(val)+"\r\\JL\\s("+totName+") Total\r\\s("+spaName+") Aggregate\r\\s("+amoName+") Amorphous\r\\s("+cGauss+") Vib Level"
	else
		KillWindow $pltName
		Display/N=$pltName/K=1 amo,spa,tot vs eV
		ModifyGraph tick=2,mirror=1,minor=1,fStyle=1
		Label left "Absorbance [a.u.]"
		Label bottom "Photon Energy [eV]"
		SetAxis bottom 1.6,3.5
		ModifyGraph lsize=1.5,rgb($spaName)=(1,26221,39321),rgb($totName)=(0,0,0)
		for(i=0;i<nGauss;i+=1)
			cGauss = StringFromList(i,gaussList)
			Wave g = $cGauss
			AppendToGraph g vs eV
			ModifyGraph lstyle($cGauss)=3,rgb($cGauss)=(21845,21845,21845),lsize($cGauss)=1.5
		endfor
		Legend/C/N=text0/J/A=RT "\\JCSet # "+num2str(val)+"\r\\JL\\s("+totName+") Total\r\\s("+spaName+") Aggregate\r\\s("+amoName+") Amorphous\r\\s("+cGauss+") Vib Level"
	endif
	SetDataFolder root:
End

Function pltSpanoSetBut(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			String panelName = "Spano_Panel"
			ControlInfo/W=$panelName fileName
			String name = S_Value		
			ControlInfo/W=$panelName setVal
			Variable val = V_Value			
			plotSpanoSet(name,val)			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarSpano1(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			String panelName = "Spano_Panel"
			ControlInfo/W=$panelName fileName
			String name = S_Value		
			ControlInfo/W=$panelName setVal
			Variable val = V_Value	
			plotSpanoSet(name,val)		
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End