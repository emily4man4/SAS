******************************************************************
*
* Program Name   : mmrmnrss.sas
* Protocol/Study : GWPharm ISE
* Type           : Analysis Reporting Dataset
* Decription     : Program to generate effect modifier analysis dataset 
*                  containing MMRM analysis of change from baseline in 
*                  NRS-spasticity score
*                  Analysis is repeated with a separate model for each 
*                  subgroup for studies GWMS0106, GWCL0403, GWSP0604 and 
*                  SAVANT 
*                  See ISE SAP Section 3.5.3
*
* Author         : Emily Foreman
* Date Created   : 08SEP2020
* Input datasets : ADSS
* Macros used    : %checklog
* Files used     : na
*
******************************************************************* 
* Change History  
*
* Changed by        : Emily Foreman
* Reason for Change : Updated subgroup categorisation variables and removed 
                      condition of at least 20 subjects in each level of a 
                      subgroup
* Date changed      : 13OCT2020
*
* Changed by        : Emily Foreman
* Reason for Change : Updated to use UN covariance structure
* Date changed      : 13NOV2020
*
* Changed by        : Emily Foreman
* Reason for Change : 1) Updated to use Part A subgroups for 0604 study
                      2) Updated to extract the p-value for the 3-way interaction term at the primary timepoint
* Date changed      : 05FEB2021
*
* Changed by        : Alice Batchelor
* Reason for Change : Added extra term needed in contrast statement for extracting the required interaction p-values from the indicator
*                     variable models coded on 02FEB. (Required adding code for two-category subgroups to also use this indicator variable
*                     method.)
* Date changed      : 26APR2021
*
* Changed by        : Alice Batchelor
* Reason for Change : Use of commas corrected in contrast statements for obtaining interaction p-values
* Date changed      : 06MAY2021
******************************************************************;

proc datasets library = work kill nolist;
run;
quit;

*Macro to test convergence of model using unstuctured covariance structure. Creates macro variable, conv, that resolves to Y if 
 the convergence criteria is met and N otherwise, in which case the compound symmetry covariance structure should be used;
%macro convergence(type = );

    %*** Check convergence to determine covariance structure ***;
    %*** ODS select statements temporarily suspend any SAS output because this MIXED procedure is only to determine convergence ***;
    ods select none;
    ods output convergencestatus = convstat;
    ods graphics off;
    proc mixed data = work.&label._adss (where = (chg ne . and trtpn ne . and &subgr. ne ' ')) method = reml;
        class trtpn usubjid avisitn (ref = "3") &subgr.;
        model chg = base trtpn avisitn &subgr. trtpn*avisitn base*avisitn &subgr.*avisitn &subgr.*trtpn &subgr.*avisitn*trtpn / ddfm = KenwardRoger; 
        repeated avisitn / type = &type. sub = usubjid;
    run;
    ods output close;
    ods select all;

	%*Use convergencestatus dataset to test convergence of the model, status = 0 and pdg = pdh = 1 indicates model has converged;
	data work.convstat1;
	  set work.convstat;
	  if status = 0 and pdg and pdh then convstat = 'Y';
	  else convstat = 'N';
	run;

%mend convergence;

*Macro to run the MMRM analysis for each study and baseline definition, i.e. repeat for secondary BL definitions for SAVANT study;
%macro mmrm_subgr(analyid  = ,
                  analydsc = ,
                  study  = ,
                  subgr  = ,
			      subgrn = ,
			      final_vis = );

  %if "&study." = "H15/02" %then %do; %let label = savant&subgrn.; %end;
  %else %do; %let label = &study.&subgrn.; %end;
  
  %*Subset data to the NRS-Spasticity Score Weekly Average for the relevant study;
  data work.&label._adss1;
    set adam.adss;
    where studyid = "&study." and paramcd = 'NRSPASAV' and basetype = "Primary";
	format &subgr.;
  run;

  %*Check that there are observation in the dataset;
  proc sql noprint;
    select *
	from work.&label._adss1;
	%let nobs = &sqlobs.;
  quit;

  %if &nobs. ne 0 %then %do;

    %********************* Create the dataset for the analysis **********************;

    %*Select observations where the relevant subgroup variable is not missing, and where, for the Savant and 0604 studies, the observation corresponds to Part B of the study;

	data work.&label._adss;
	  set work.&label._adss1;
	  where &subgr. ne ' ' and lowcase(aphase) = 'double blind treatment';
	run;

	%******************** Calculate the data for the analysis reporting dataset ************************;
    
    %*Calculate the small n counts overall and for each week, in each treatment group, for each subgroup. The MMRM model excludes subjects with missing covariates,
      therefore need to also exclude these from the small n count;

	proc freq data = work.&label._adss (where = (chg ne . and trtpn ne . and &subgr. ne ' '));
	  tables trtpn*avisitn*&subgr. / out = work.&label._counts (drop = percent);
	run;

	data work.&label._counts1 (keep = trtn &subgr. timeptw n);
	  length &subgr. timeptw $ 200;
      set work.&label._counts (where = (avisitn = &final_vis.)
		                       rename = (trtpn = trtn
								         count = n));
    timeptw = "Week &final_vis.";
	run;
      
    %*Call the convergence macro as defined above to test convergence using unstructured covariance structure;
    %convergence(type = un); 

	%*Create macro variable to indicate model convergence using unstructured covariance structure;
    data _null_;
	  set work.convstat1;
	  call symputx("conv", convstat);
	run;

	%put &conv.;
      
    %*If model converges then covariance structure is unstructured
     else if model does not converge then check whether model converges using type = cs;
    %if "&conv." = "Y" %then %do;

      %let type = un; 
      
      ods output lsmeans = work.&label._lsmeans diffs = work.&label._diffs solutionF = work.&label._soln;
      proc mixed data = work.&label._adss (where = (chg ne . and trtpn ne . and &subgr. ne ' ')) method = reml;
        class trtpn avisitn (ref = "1") usubjid &subgr.;
        model chg = base trtpn avisitn &subgr. trtpn*avisitn base*avisitn &subgr.*avisitn &subgr.*trtpn &subgr.*avisitn*trtpn / s DDFM = KenwardRoger;
        repeated avisitn / type = un subject = usubjid rcorr;  
        lsmeans trtpn*avisitn*&subgr. / cl pdiff;
      run;

      data work.&label._lsmeans;
		length covstr $ 200;
		set work.&label._lsmeans;
		covstr = 'un'; 
      run;

      %*Create macro variable to indicate whether there was enough data for model to converge;
      %let converge = Y;
      
    %end;
    %else %do;

	  %convergence(type = cs);

      %*Create macro variable to indicate model convergence using compound symmetry covariance structure;
      data _null_;
	    set work.convstat1;
	    call symputx("conv", convstat);
	  run;

      %put &conv.;

      %*If conv = Y then type = cs else test convergence using the ar(1) covariance strucutre;
      %if "&conv." = "Y" %then %do;

        %let type = cs;
      
        ods output lsmeans = work.&label._lsmeans diffs = work.&label._diffs solutionF = work.&label._soln;
        proc mixed data = work.&label._adss (where = (chg ne . and trtpn ne . and &subgr. ne ' ')) method = reml;
          class trtpn avisitn (ref = "1") usubjid &subgr.;
          model chg = base trtpn avisitn &subgr. trtpn*avisitn base*avisitn &subgr.*avisitn &subgr.*trtpn &subgr.*avisitn*trtpn / s DDFM = KenwardRoger;
          repeated avisitn / type = cs subject = usubjid rcorr;
          lsmeans trtpn*avisitn*&subgr. / cl pdiff;
        run; 

		data work.&label._lsmeans;
		  length covstr $ 200;
		  set work.&label._lsmeans;
		  covstr = 'cs';
		run;

		%*Create macro variable to indicate whether there was enough data for model to converge;
		%let converge = Y;

      %end;
      %else %do;

	    %convergence(type = %str(ar(1)));
		 
		%*Create macro variable to indicate model convergence using ar(1) covariance structure;
        data _null_;
	      set work.convstat1;
	      call symputx("conv", convstat);
	    run;
		  
		%put &conv.;

		%*If conv = Y then type = ar(1) else set converge macro variable to N to indicate insufficient data for model to converge;
		%if "&conv." = "Y" %then %do;

          %let type = %str(ar(1));

          ods output lsmeans = work.&label._lsmeans diffs = work.&label._diffs solutionF = work.&label._soln;
          proc mixed data = work.&label._adss (where = (chg ne . and trtpn ne . and &subgr. ne ' ')) method = reml;
            class trtpn avisitn (ref = "1") usubjid &subgr.;
            model chg = base trtpn avisitn &subgr. trtpn*avisitn base*avisitn &subgr.*avisitn &subgr.*trtpn &subgr.*avisitn*trtpn / s DDFM = KenwardRoger;
            repeated avisitn / type = ar(1) subject = usubjid rcorr;
            lsmeans trtpn*avisitn*&subgr. / cl pdiff;
          run; 
		    
		  data work.&label._lsmeans;
		    length covstr $ 200;
		    set work.&label._lsmeans;
		    covstr = 'ar(1)';
		  run;
		    
		  %*Create macro variable to indicate whether there was enough data for model to converge;
		  %let converge = Y;

		%end;
		%else %do;

		  %let converge = N;

		%end;

      %end;
      
    %end;

	%if "&converge." = "Y" %then %do;

	  proc sql noprint;
		create table work.&label._catcount as
		select distinct &subgr.
        from work.&label._adss
	    where &subgr. ne "";
		select count(*) as N
		into: cat_count
		from work.&label._catcount;
	  quit;

	  %*If the current subgroup has >2 categories (agegr2 and conspan1) then create indicator variables to fit the model and use a contrast statement to 
		extract the 3-way interaction p-value at the primary timepoint;
      %*AB 26APR2021: Fit model using indicator variables for 2-category subgroups too (extra term required to find p-values);

      %if "&subgr." ne "agegr2" and "&subgr." ne "conspan1" %then %do;                          
	    data work.&label._adss2;
          set work.&label._adss;    

		          if avisitn = 2 then vis2 = 1;
		          else vis2 = 0;
		          if avisitn = 3 then vis3 = 1;
		          else vis3 = 0;
		          if avisitn = 4 then vis4 = 1;
		          else vis4 = 0;
		          if avisitn = 5 then vis5 = 1;
		          else vis5 = 0;
		          if avisitn = 6 then vis6 = 1;
		          else vis6 = 0;
		          %if "&study." ne "GWMS0106" %then %do;
		            if avisitn = 7 then vis7 = 1;
		            else vis7 = 0;
		            if avisitn = 8 then vis8 = 1;
		            else vis8 = 0;
		            if avisitn = 9 then vis9 = 1;
		            else vis9 = 0;
		            if avisitn = 10 then vis10 = 1;
		            else vis10 = 0;
		            if avisitn = 11 then vis11 = 1;
		            else vis11 = 0;
		            if avisitn = 12 then vis12 = 1;
		            else vis12 = 0;
		          %end;
		          %if "&study." = "GWCL0403" %then %do;
		            if avisitn = 13 then vis13 = 1;
		            else vis13 = 0;
		            if avisitn = 14 then vis14 = 1;
		            else vis14 = 0;
		          %end;

                  %if "&subgr." = "agegr1" %then %do;
                    if &subgr.='>=50' then n1=1;
                  %end;
                  %else %if "&subgr." = "ashbl1" %then %do;
                    if &subgr.='>=2.3' then n1=1;
                  %end;
                  %else %if "&subgr." = "baclfus1" or "&subgr." = "bndpuse1" or "&subgr." = "tznduse1" %then %do;                          
                    if &subgr.='Not Currently Taking' then n1=1;
                  %end;
                  %else %if "&subgr." = "canusefl" %then %do;
                    if &subgr.='Y' then n1=1;                               
                  %end;
                  %else %if "&subgr." = "consusb1" or "&subgr." = "consusn1" %then %do;                          
                    if &subgr.='Yes' then n1=1;                                    
                  %end;
                  %else %if "&subgr." = "sex" %then %do;
                    if &subgr.='M' then n1=1;                             
                  %end;
                  %else %if "&subgr." = "nrsabl1" or "&subgr." = "nrsbl1" %then %do;                          
                    if &subgr.='>=4' then n1=1;
                  %end;
                  %else %if "&subgr." = "nrsabl2" or "&subgr." = "nrsbl2" %then %do;                        
                    if &subgr.='>=5' then n1=1;
                  %end;
                  %else %if "&subgr." = "nrsabl3" or "&subgr." = "nrsbl3" or "&subgr." = "edsslbl" %then %do;                         
                    if &subgr.='>=6' then n1=1;
                  %end;
	              %else %if "&subgr." = "nrsabl4" or "&subgr." = "nrsbl4" %then %do;                         
                    if &subgr.='>=7' then n1=1;
                  %end;
                  %else %if "&subgr." = "nrsabl5" or "&subgr." = "nrsbl5" %then %do;                        
                    if &subgr.='>=8' then n1=1;
                  %end;
                  %else %if "&subgr." = "masabl1" or "&subgr." = "masbl1" %then %do;                         
                    if &subgr.='>=1' then n1=1;
                  %end;
                  %else %if "&subgr." = "mssdur1" %then %do;
                    if &subgr.='>=7 years' then n1=1;
                  %end;

                  else n1=0;
            run;

		%if "&study." = "GWMS0106" %then %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 vis2 vis3 vis4 vis5 vis6
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6
			  	        n1*trtpn 
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
            contrast "Treatment*Week 6*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis6 1; %* AB 06MAY2021: Removed comma between terms;
          run;

		%end;
		%else %if "&study." = "GWCL0403" %then %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 vis2 vis3 vis4 vis5 vis6 vis7 vis8 vis9 vis10 vis11 vis12 vis13 vis14
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6 trtpn*vis7 trtpn*vis8 
                        trtpn*vis9 trtpn*vis10 trtpn*vis11 trtpn*vis12 trtpn*vis13 trtpn*vis14
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6 base*vis7 base*vis8 
                        base*vis9 base*vis10 base*vis11 base*vis12 base*vis13 base*vis14
			  	        n1*trtpn 
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6 n1*vis7 n1*vis8 n1*vis9 n1*vis10 n1*vis11 n1*vis12 n1*vis13 n1*vis14
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 n1*trtpn*vis7 n1*trtpn*vis8 
                        n1*trtpn*vis9 n1*trtpn*vis10 n1*trtpn*vis11 n1*trtpn*vis12 n1*trtpn*vis13 n1*trtpn*vis14 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
            contrast "Treatment*Week 14*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis14 1; %* AB 06MAY2021: Removed comma between terms;
          run;

		%end;
		%else %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 vis2 vis3 vis4 vis5 vis6 vis7 vis8 vis9 vis10 vis11 vis12
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6 trtpn*vis7 trtpn*vis8 trtpn*vis9 trtpn*vis10 trtpn*vis11 trtpn*vis12
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6 base*vis7 base*vis8 base*vis9 base*vis10 base*vis11 base*vis12
			  	        n1*trtpn 
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6 n1*vis7 n1*vis8 n1*vis9 n1*vis10 n1*vis11 n1*vis12
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 n1*trtpn*vis7 
                        n1*trtpn*vis8 n1*trtpn*vis9 n1*trtpn*vis10 n1*trtpn*vis11 n1*trtpn*vis12 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
            contrast "Treatment*Week 12*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis12 1; %* AB 06MAY2021: Removed comma between terms;
          run;

		%end;

		data _null_;
		  set work.&label._cont;
		  call symputx("pvalue", probf);
        run;

      %end;

	  %*------ AB 26APR2021: End of new code ----------;

	  %else %if "&subgr." = "agegr2" %then %do; %* AB 26APR2021: Updated if to else if;

	    data work.&label._adss2;
		  set work.&label._adss;
		  if avisitn = 2 then vis2 = 1;
		  else vis2 = 0;
		  if avisitn = 3 then vis3 = 1;
		  else vis3 = 0;
		  if avisitn = 4 then vis4 = 1;
		  else vis4 = 0;
		  if avisitn = 5 then vis5 = 1;
		  else vis5 = 0;
		  if avisitn = 6 then vis6 = 1;
		  else vis6 = 0;
		  %if "&study." ne "GWMS0106" %then %do;
		    if avisitn = 7 then vis7 = 1;
		    else vis7 = 0;
		    if avisitn = 8 then vis8 = 1;
		    else vis8 = 0;
		    if avisitn = 9 then vis9 = 1;
		    else vis9 = 0;
		    if avisitn = 10 then vis10 = 1;
		    else vis10 = 0;
		    if avisitn = 11 then vis11 = 1;
		    else vis11 = 0;
		    if avisitn = 12 then vis12 = 1;
		    else vis12 = 0;
		  %end;
		  %if "&study." = "GWCL0403" %then %do;
		    if avisitn = 13 then vis13 = 1;
		    else vis13 = 0;
		    if avisitn = 14 then vis14 = 1;
		    else vis14 = 0;
		  %end;
		  if agegr2 = ">40-50" then n1 = 1;
          else n1 = 0;
	      if agegr2 = ">50-60" then n2 = 1;
          else n2 = 0;
		  if agegr2 = ">60" then n3 = 1;
          else n3 = 0;
		run;

		%if "&study." = "GWMS0106" %then %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 n2 n3 vis2 vis3 vis4 vis5 vis6
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6
			  	        n1*trtpn n2*trtpn n3*trtpn
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6
			  		    n2*vis2 n2*vis3 n2*vis4 n2*vis5 n2*vis6
			  		    n3*vis2 n3*vis3 n3*vis4 n3*vis5 n3*vis6
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6
			  		    n2*trtpn*vis2 n2*trtpn*vis3 n2*trtpn*vis4 n2*trtpn*vis5 n2*trtpn*vis6
			  		    n3*trtpn*vis2 n3*trtpn*vis3 n3*trtpn*vis4 n3*trtpn*vis5 n3*trtpn*vis6 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
			%* AB 26APR2021: Extra required 2 factor interaction terms added to contrast statement;
			%* AB 06MAY2021: Commas only between each subgroup category;
            contrast "Treatment*Week 6*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis6 1, n2*trtpn 1 n2*trtpn*vis6 1, n3*trtpn 1 n3*trtpn*vis6 1;
          run;

		%end;
		%else %if "&study." = "GWCL0403" %then %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 n2 n3 vis2 vis3 vis4 vis5 vis6 vis7 vis8 vis9 vis10 vis11 vis12 vis13 vis14
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6 trtpn*vis7 trtpn*vis8 
                        trtpn*vis9 trtpn*vis10 trtpn*vis11 trtpn*vis12 trtpn*vis13 trtpn*vis14
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6 base*vis7 base*vis8 
                        base*vis9 base*vis10 base*vis11 base*vis12 base*vis13 base*vis14
			  	        n1*trtpn n2*trtpn n3*trtpn
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6 n1*vis7 n1*vis8 n1*vis9 n1*vis10 n1*vis11 n1*vis12 n1*vis13 n1*vis14
			  		    n2*vis2 n2*vis3 n2*vis4 n2*vis5 n2*vis6 n2*vis7 n2*vis8 n2*vis9 n2*vis10 n2*vis11 n2*vis12 n2*vis13 n2*vis14
			  		    n3*vis2 n3*vis3 n3*vis4 n3*vis5 n3*vis6 n3*vis7 n3*vis8 n3*vis9 n3*vis10 n3*vis11 n3*vis12 n3*vis13 n3*vis14
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 n1*trtpn*vis7 n1*trtpn*vis8 
                        n1*trtpn*vis9 n1*trtpn*vis10 n1*trtpn*vis11 n1*trtpn*vis12 n1*trtpn*vis13 n1*trtpn*vis14
			  		    n2*trtpn*vis2 n2*trtpn*vis3 n2*trtpn*vis4 n2*trtpn*vis5 n2*trtpn*vis6 n2*trtpn*vis7 n2*trtpn*vis8 
                        n2*trtpn*vis9 n2*trtpn*vis10 n2*trtpn*vis11 n2*trtpn*vis12 n2*trtpn*vis13 n2*trtpn*vis14
			  		    n3*trtpn*vis2 n3*trtpn*vis3 n3*trtpn*vis4 n3*trtpn*vis5 n3*trtpn*vis6 n3*trtpn*vis7 n3*trtpn*vis8 
                        n3*trtpn*vis9 n3*trtpn*vis10 n3*trtpn*vis11 n3*trtpn*vis12 n3*trtpn*vis13 n3*trtpn*vis14 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
			%* AB 26APR2021: Extra required 2 factor interaction terms added to contrast statement;
			%* AB 06MAY2021: Commas only between each subgroup category;
            contrast "Treatment*Week 14*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis14 1, n2*trtpn 1 n2*trtpn*vis14 1, n3*trtpn 1 n3*trtpn*vis14 1;
          run;

		%end;
		%else %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 n2 n3 vis2 vis3 vis4 vis5 vis6 vis7 vis8 vis9 vis10 vis11 vis12
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6 trtpn*vis7 trtpn*vis8 trtpn*vis9 trtpn*vis10 trtpn*vis11 trtpn*vis12
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6 base*vis7 base*vis8 base*vis9 base*vis10 base*vis11 base*vis12
			  	        n1*trtpn n2*trtpn n3*trtpn
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6 n1*vis7 n1*vis8 n1*vis9 n1*vis10 n1*vis11 n1*vis12
			  		    n2*vis2 n2*vis3 n2*vis4 n2*vis5 n2*vis6 n2*vis7 n2*vis8 n2*vis9 n2*vis10 n2*vis11 n2*vis12
			  		    n3*vis2 n3*vis3 n3*vis4 n3*vis5 n3*vis6 n3*vis7 n3*vis8 n3*vis9 n3*vis10 n3*vis11 n3*vis12 
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 n1*trtpn*vis7 
                        n1*trtpn*vis8 n1*trtpn*vis9 n1*trtpn*vis10 n1*trtpn*vis11 n1*trtpn*vis12
			  		    n2*trtpn*vis2 n2*trtpn*vis3 n2*trtpn*vis4 n2*trtpn*vis5 n2*trtpn*vis6 n2*trtpn*vis7 
                        n2*trtpn*vis8 n2*trtpn*vis9 n2*trtpn*vis10 n2*trtpn*vis11 n2*trtpn*vis12
			  		    n3*trtpn*vis2 n3*trtpn*vis3 n3*trtpn*vis4 n3*trtpn*vis5 n3*trtpn*vis6 n3*trtpn*vis7 
                        n3*trtpn*vis8 n3*trtpn*vis9 n3*trtpn*vis10 n3*trtpn*vis11 n3*trtpn*vis12 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
			%* AB 26APR2021: Extra required 2 factor interaction terms added to contrast statement;
			%* AB 06MAY2021: Commas only between each subgroup category;
            contrast "Treatment*Week 12*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis12 1, n2*trtpn 1 n2*trtpn*vis12 1, n3*trtpn 1 n3*trtpn*vis12 1;
          run;

		%end;

	    data _null_;
		  set work.&label._cont;
		  call symputx("pvalue", probf);
        run;

	  %end;
	  %else %if "&subgr." = "conspan1" %then %do;

		data work.&label._adss2;
		  set work.&label._adss;
		  if avisitn = 2 then vis2 = 1;
		  else vis2 = 0;
		  if avisitn = 3 then vis3 = 1;
		  else vis3 = 0;
		  if avisitn = 4 then vis4 = 1;
		  else vis4 = 0;
		  if avisitn = 5 then vis5 = 1;
		  else vis5 = 0;
		  if avisitn = 6 then vis6 = 1;
		  else vis6 = 0;
		  %if "&study." ne "GWMS0106" %then %do;
		    if avisitn = 7 then vis7 = 1;
		    else vis7 = 0;
		    if avisitn = 8 then vis8 = 1;
		    else vis8 = 0;
		    if avisitn = 9 then vis9 = 1;
		    else vis9 = 0;
		    if avisitn = 10 then vis10 = 1;
		    else vis10 = 0;
		    if avisitn = 11 then vis11 = 1;
		    else vis11 = 0;
		    if avisitn = 12 then vis12 = 1;
		    else vis12 = 0;
		  %end;
		  %if "&study." = "GWCL0403" %then %do;
		    if avisitn = 13 then vis13 = 1;
		    else vis13 = 0;
		    if avisitn = 14 then vis14 = 1;
		    else vis14 = 0;
		  %end;
		  if conspan1 = "1" then n1 = 1;
          else n1 = 0;
		  if conspan1 = "2+" then n2 = 1;
          else n2 = 0;
		run;

		%if "&study." = "GWMS0106" %then %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 n2 vis2 vis3 vis4 vis5 vis6
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6
			  	        n1*trtpn n2*trtpn
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6
			  		    n2*vis2 n2*vis3 n2*vis4 n2*vis5 n2*vis6
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6
			  		    n2*trtpn*vis2 n2*trtpn*vis3 n2*trtpn*vis4 n2*trtpn*vis5 n2*trtpn*vis6/ s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
			%* AB 26APR2021: Extra required 2 factor interaction terms added to contrast statement;
			%* AB 06MAY2021: Commas only between each subgroup category;
            contrast "Treatment*Week 6*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis6 1, n2*trtpn 1 n2*trtpn*vis6 1;
          run;

		%end;
		%else %if "&study." = "GWCL0403" %then %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 n2 vis2 vis3 vis4 vis5 vis6 vis7 vis8 vis9 vis10 vis11 vis12 vis13 vis14
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6 trtpn*vis7 trtpn*vis8 
                        trtpn*vis9 trtpn*vis10 trtpn*vis11 trtpn*vis12 trtpn*vis13 trtpn*vis14
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6 base*vis7 base*vis8 
                        base*vis9 base*vis10 base*vis11 base*vis12 base*vis13 base*vis14
			  	        n1*trtpn n2*trtpn
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6 n1*vis7 n1*vis8 n1*vis9 n1*vis10 n1*vis11 n1*vis12 n1*vis13 n1*vis14
			  		    n2*vis2 n2*vis3 n2*vis4 n2*vis5 n2*vis6 n2*vis7 n2*vis8 n2*vis9 n2*vis10 n2*vis11 n2*vis12 n2*vis13 n2*vis14
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 n1*trtpn*vis7 n1*trtpn*vis8 
                        n1*trtpn*vis9 n1*trtpn*vis10 n1*trtpn*vis11 n1*trtpn*vis12 n1*trtpn*vis13 n1*trtpn*vis14
			  		    n2*trtpn*vis2 n2*trtpn*vis3 n2*trtpn*vis4 n2*trtpn*vis5 n2*trtpn*vis6 n2*trtpn*vis7 n2*trtpn*vis8 
                        n2*trtpn*vis9 n2*trtpn*vis10 n2*trtpn*vis11 n2*trtpn*vis12 n2*trtpn*vis13 n2*trtpn*vis14 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
			%* AB 26APR2021: Extra required 2 factor interaction terms added to contrast statement;
			%* AB 06MAY2021: Commas only between each subgroup category;
            contrast "Treatment*Week 14*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis14 1, n2*trtpn 1 n2*trtpn*vis14 1;
          run;

		%end;
		%else %do;
		    
          ods output contrasts = work.&label._cont;
          proc mixed data = work.&label._adss2 (where = (chg ne . and trtpn ne . and &subgr. ne '')) method = reml;
            class avisitn usubjid;
            model chg = base trtpn n1 n2 vis2 vis3 vis4 vis5 vis6 vis7 vis8 vis9 vis10 vis11 vis12
			            trtpn*vis2 trtpn*vis3 trtpn*vis4 trtpn*vis5 trtpn*vis6 trtpn*vis7 trtpn*vis8 trtpn*vis9 trtpn*vis10 trtpn*vis11 trtpn*vis12
                        base*vis2 base*vis3 base*vis4 base*vis5 base*vis6 base*vis7 base*vis8 base*vis9 base*vis10 base*vis11 base*vis12
			  	        n1*trtpn n2*trtpn
			  		    n1*vis2 n1*vis3 n1*vis4 n1*vis5 n1*vis6 n1*vis7 n1*vis8 n1*vis9 n1*vis10 n1*vis11 n1*vis12
			  		    n2*vis2 n2*vis3 n2*vis4 n2*vis5 n2*vis6 n2*vis7 n2*vis8 n2*vis9 n2*vis10 n2*vis11 n2*vis12
			  		    n1*trtpn*vis2 n1*trtpn*vis3 n1*trtpn*vis4 n1*trtpn*vis5 n1*trtpn*vis6 n1*trtpn*vis7 
                        n1*trtpn*vis8 n1*trtpn*vis9 n1*trtpn*vis10 n1*trtpn*vis11 n1*trtpn*vis12
			  		    n2*trtpn*vis2 n2*trtpn*vis3 n2*trtpn*vis4 n2*trtpn*vis5 n2*trtpn*vis6 n2*trtpn*vis7 
                        n2*trtpn*vis8 n2*trtpn*vis9 n2*trtpn*vis10 n2*trtpn*vis11 n2*trtpn*vis12 / s DDFM = KenwardRoger;
            repeated avisitn / type = &type. subject = usubjid;
			%* AB 26APR2021: Extra required 2 factor interaction terms added to contrast statement;
			%* AB 06MAY2021: Commas only between each subgroup category;
            contrast "Treatment*Week 12*Subgroup Interaction" n1*trtpn 1 n1*trtpn*vis12 1, n2*trtpn 1 n2*trtpn*vis12 1;
          run;

		%end;

	    data _null_;
		  set work.&label._cont;
		  call symputx("pvalue", probf);
        run;

	  %end;

      %* AB 26APR2021: Removed code using fixed effect solutions for 2-category subgroups p-values since no longer needed,
	                   added check for missing p-value here;

      %if &pvalue. = %then %do;
	    %let pvalue = .;
	  %end;
      
	  %*Format the least squares estimates and then merge them onto the small n counts;
      data work.&label._lsmeans1 (keep = &subgr. timeptw timeptwn trt trtn lsmean lsmeanse lowerci upperci covstr);
        length &subgr. timeptw $ 200;
        set work.&label._lsmeans (rename = (trtpn = trtn
                                            estimate = lsmean
	      									stderr = lsmeanse
	      									lower = lowerci
	      									upper = upperci));
		where avisitn = &final_vis.;
	    timeptw = "Week &final_vis.";
		timeptwn = &final_vis.;
	    if trtn = 1 then trt = 'Nabiximols';
	    else if trtn = 2 then trt = 'Placebo';
      run;
        
      proc sql;
        create table work.&label._lsmeans_merge as
	    select one.*, two.n
	    from work.&label._lsmeans1 as one left join work.&label._counts1 as two
	    on one.trtn = two.trtn and one.&subgr. = two.&subgr.;
      quit;
        
	  %*Format the least squares mean difference estimates and then merge them onto the least squares dataset;
      data work.&label._diffs1 (keep = trtn &subgr. lsmdif lsmdlci lsmduci pval);
	    length &subgr. $ 200;
        set work.&label._diffs (rename = (trtpn = trtn
	                                      _trtpn = _trtn
	                                      estimate = lsmdif
	      								lower = lsmdlci
	      								upper = lsmduci
	      								probt = pval));
        where avisitn = &final_vis. and trtn = 1 and _trtn = 2 and avisitn = _avisitn and &subgr. = _&subgr.;
      run;
	    
      proc sql;
        create table work.&label._diffs_merge (rename = (&subgr. = category)) as
        select one.*, two.lsmdif, two.lsmdlci, two.lsmduci, two.pval
	    from work.&label._lsmeans_merge as one left join work.&label._diffs1 as two
	    on one.trtn = two.trtn and one.&subgr. = two.&subgr.;
      quit;
	    
	  data work.&label._tests1;
	    pval = &pvalue.;
	  run;
	    
	  %*Format the merged dataset, including an order variable for the categories;
	  data work.&label._merge;
	    length category $ 200;
	    set work.&label._tests1 (in = A)
        work.&label._diffs_merge (in = B);
        analyid = "&analyid.";
		analydsc = "&analydsc.";
	    studyid = "&study.";
		subgroup = "&subgr.";
	    if A then do;
          category = "Treatment*Week &final_vis.*&subgr. Interaction";
		  categoryn = 2;
	      timeptw = "Week &final_vis.";
		  timeptwn = &final_vis.;
	    end;
		else if B then do;
          category = category;
		  categoryn = 1;
		end;
	  run;

	%end;

	%*If there was not enough data for the model to converge using UN, CS or AR(1) then create dataset containing the counts only;
	%else %if "&converge." = "N" %then %do;

	  proc freq data = work.&label._adss (where = (chg ne . and trtpn ne . and &subgr. ne ' '));
	    tables trtpn*avisitn*&subgr. / out = work.&label._counts (drop = percent);
	  run;
	    
	  %*Create a shell, with the distinct levels of the subgroup variable that can be observed, repeated for each treatment;
      proc sort data = work.adss_merge (where = (&subgr. ne ' ')
                                        keep = &subgr.)
                 out = work.&label._shell nodupkey;
        by &subgr.;
      run;
        
      data work.&label._shell;
        set work.&label._shell;
	    do trtpn = 1 to 2;
	      output;
	    end;
      run;
	    
	  %*Merge the subjects with the shell to ensure each level of the subgroup is accounted for;
      proc sql;
        create table work.&label._counts_merge as
	    select one.*, two.count, two.avisitn
	    from work.&label._shell as one left join work.&label._counts (where = (avisitn = &final_vis.)) as two
	    on one.trtpn = two.trtpn and one.&subgr. = two.&subgr.;
      quit;
	    
	  data work.&label._merge (rename = (&subgr. = category trtpn = trtn)
                               drop = avisitn);
	    length &subgr. timeptw covstr $ 200;
		set work.&label._counts_merge (where = (avisitn = &final_vis. or avisitn = .)
                                       rename = (count = n));
		if n = . then n = 0;
		analyid = "&analyid.";
		analydsc = "&analydsc.";
		studyid = "&study.";
		subgroup = "&subgr.";
		categoryn = 1;
		if trtpn = 1 then trt = 'Nabiximols';
		else trt = 'Placebo';
		timeptw = "Week &final_vis.";
		timeptwn = &final_vis.;
		covstr = '';
	    lsmean = .;
	    lsmeanse = .;
	    lowerci = .;
	    upperci = .;
        lsmdif = .;
	    lsmdlci = .;
	    lsmduci = .;
	    pval = .;
	  run;

	%end;

  %end;

  %*Create empty dataset if there were no observations in the dataset;
  %else %do;

    data work.&label._merge;
	  length category trt timeptw covstr $ 200;
      analyid = "&analyid.";
	  analydsc = "&analydsc.";
      studyid = "&study.";
      subgroup = "&subgr.";
	  category = ' ';
	  categoryn = .;
      trt = '';
	  trtn = .;
      timeptw = '';
	  timeptwn = .;
	  covstr = '';
	  n = .;
	  lsmean = .;
	  lsmeanse = .;
	  lowerci = .;
	  upperci = .;
	  lsmdif = .;
	  lsmdlci = .;
	  lsmduci = .;
	  pval = .;
    run;

  %end;

%mend mmrm_subgr;

%*Macro_run is a macro used to call the mmrm_subgr analysis for each study, and repeat over all subgroups;

*********
Subgroup Variables:

- Age Groups: AGEGR1, AGEGR2
- Sex: SEX
- Baseline EDSS: EDSSLBL (excluding GWMS0106)
- Prior Cannabis use: CANUSEFL (excluding Savant)
- Concomitant anti-spasticity medication use (narrow definition) : CONSUSN1
- Concomitant anti-spasticity medication use (broader definition): CONSUSB1
- Number of concomitant anti-spasticity medications: CONSPAN1
- Baclofen use: BACLFUS1
- Tizanidine use: TZNDUSE1
- Benzodiazepine derivative use: BNDPUSE1
- Severity of illness at BL based on initial trial inclusion BL NRS-spasticity score (<4, >=4): NRSBL1 (excluding GWSP0604)
- Severity of illness at BL based on initial trial inclusion BL NRS-spasticity score (<5, >=5): NRSBL2 (excluding GWSP0604)
- Severity of illness at BL based on initial trial inclusion BL NRS-spasticity score (<6, >=6): NRSBL3
- Severity of illness at BL based on initial trial inclusion BL NRS-spasticity score (<7, >=7): NRSBL4
- Severity of illness at BL based on initial trial inclusion BL NRS-spasticity score (<8, >=8): NRSBL5 (GWSP0604 only)
- Severity of illness at BL based on duration of spasticity symptoms (<7yrs, >=7 yrs): MSSDUR1 (excluding GWMS0106)
- Baseline Ashworth Score: ASHBL1 (GWMS0106 only)
- Modified Ashworth Score: MASBL1 (excluding GWMS0106)
*;

%macro macro_run(anlid = ,
                 anlydsc = ,
                 stud = ,
				 fvis = );

  %*Define macro variables for each subgroup;
  %let subgr1 = agegr1;
  %let subgr2 = agegr2;
  %let subgr3 = sex; 
  %let subgr4 = edsslbl;
  %let subgr5 = canusefl;
  %let subgr6 = consusn1;
  %let subgr7 = consusb1;
  %let subgr8 = conspan1;
  %let subgr9 = baclfus1;
  %let subgr10 = tznduse1;
  %let subgr11 = bndpuse1;
  %let subgr12 = nrsbl1;
  %let subgr13 = nrsbl2; 
  %let subgr14 = nrsbl3; 
  %let subgr15 = nrsbl4; 
  %let subgr16 = nrsbl5; 
  %let subgr17 = mssdur1; 
  %let subgr18 = ashbl1;  
  %let subgr19 = masbl1; 

  %*Use Part A categorizations for the effect modifier subgroups for efficacy endpoints in 0604 study;
  %if "&stud." = "GWSP0604" %then %do;

    %let subgr12 = nrsabl1;
    %let subgr13 = nrsabl2; 
    %let subgr14 = nrsabl3; 
    %let subgr15 = nrsabl4; 
    %let subgr16 = nrsabl5;  
    %let subgr19 = masabl1; 

  %end;  

  %do i = 1 %to 19;
  
    %*Data not collected for BL EDSS, duration of spasticity symptoms at BL, modified ashworth score or BL NRS Score Group 5 for 0106 study therefore do not run analysis for corresponding subgroup variables;
    %if "&stud." = "GWMS0106" %then %do;
      %if &i. ne 4 and &i. ne 16 and &i. ne 17 and &i. ne 19 %then %do; 
        %mmrm_subgr(analyid = &anlid., analydsc = &anlydsc., study = &stud., subgr = &&subgr&i.., subgrn = &i., final_vis = &fvis.)
      %end;
    %end;

	%*Data for baseline ashworth score not collected for GWCL0403, GWSP0604 and savant studies therefore do not run analysis for corresponding subgroup variable;
	%else %if &i. ne 18  %then %do;
	  %if "&stud." = "H15/02" %then %do;
	    %*Data for prior cannabis use and BL NRS Score Group 5 not collected for savant study therefore do not run analysis for corresponding subgroup variable;
	    %if &i. ne 5 and &i. ne 16 %then %do;
          %mmrm_subgr(analyid = &anlid., analydsc = &anlydsc., study = &stud., subgr = &&subgr&i.., subgrn = &i., final_vis = &fvis.)
	    %end;
	  %end;

	  %*Data for BL NRS Score Group 5 not collected for 0403 study therefore do not run analysis for corresponding subgroup variable;
      %else %if "&stud." = "GWCL0403" %then %do;
	    %if &i. ne 16 %then %do;
	      %mmrm_subgr(analyid = &anlid., analydsc = &anlydsc., study = &stud., subgr = &&subgr&i.., subgrn = &i., final_vis = &fvis.)
	    %end;
      %end;

	  %*Data for BL NRS Score Groups 1 and 2 not collected for 0604 study therefore do not run analysis for corresponding subgroup variable;
      %else %if "&stud." = "GWSP0604" %then %do;
	    %if &i. ne 12 and &i. ne 13 %then %do;
	      %mmrm_subgr(analyid = &anlid., analydsc = &anlydsc., study = &stud., subgr = &&subgr&i.., subgrn = &i., final_vis = &fvis.)
	    %end;
      %end;
	%end;
  
  %end;

  %*Merge subgroup analysis datasets into overall datasets for the corresponding study;

  %if "&stud." = "GWMS0106" %then %do;

    data work.&stud._merge;
	  length subgroup category trt $ 200;
      set %do i = 1 %to 19;
              %if &i. ne 4 and &i. ne 16 and &i. ne 17 and &i. ne 19 %then %do;
	  		    work.&stud.&i._merge
	  		  %end;
	  	  %end;
	  	  ;
    run;

  %end;
  %else %if "&stud." = "H15/02" %then %do; 

    data work.savant_merge;
	  length subgroup category trt $ 200;
      set %do i = 1 %to 19;
              %if &i. ne 5 and &i. ne 16 and &i. ne 18 %then %do;
	  		    work.savant&i._merge
	  		  %end;
	  	  %end;
	  	  ;
    run;

  %end;
  %else %if "&stud." = "GWCL0403" %then %do;

    data work.&stud._merge;
	  length subgroup category trt $ 200;
	  set %do i = 1 %to 19;
	        %if &i. ne 16 and &i. ne 18 %then %do;
	          work.&stud.&i._merge
			%end;
		  %end;
		  ;
	run;

  %end;
  %else %if "&stud." = "GWSP0604" %then %do;

    data work.&stud._merge;
	  length subgroup category trt $ 200;
	  set %do i = 1 %to 19;
	        %if &i. ne 12 and &i. ne 13 and &i. ne 18 %then %do;
	          work.&stud.&i._merge
			%end;
		  %end;
		  ;
	run;

  %end;

%mend macro_run;

%macro_run(anlid = mmrmemnrs01, anlydsc = %str(MMRM Effect Modification of CFB in NRS-Spasticity by Subgroup During the Last Week of Treatment), stud = GWMS0106, fvis = 6)
%macro_run(anlid = mmrmemnrs02, anlydsc = %str(MMRM Effect Modification of CFB in NRS-Spasticity by Subgroup During the Last Week of Treatment), stud = GWCL0403, fvis = 14)
%macro_run(anlid = mmrmemnrs03, anlydsc = %str(MMRM Effect Modification of CFB in NRS-Spasticity by Subgroup During the Last Week of Treatment), stud = GWSP0604, fvis = 12)
%macro_run(anlid = mmrmemnrs04, anlydsc = %str(MMRM Effect Modification of CFB in NRS-Spasticity by Subgroup During the Last Week of Treatment), stud = H15/02,   fvis = 12)

%*Merge the data for each study into one output dataset;
data work.mmrmnrss;
attrib
    analyid   label = "Analysis ID"
    analydsc  label = "Analysis Description"
	pop       label = "Analysis Set"
	studyid   label = "Study Identifier"
	subgroup  label = "Subgroup"
	category  label = "Category"
	trt       label = "Planned Treatment"
	trtn      label = "Planned Treatment (N)"
	timeptw   label = "Analysis Timepoint"
	timeptwn  label = "Analysis Timepoint (N)"
	covstr    label = "Covariance Structure"
	n         label = "n"
	lsmean    label = "LS Mean"
	lsmeanse  label = "SE of LS Mean"
	lowerci   label = "LS Mean Low CI"
	upperci   label = "LS Mean Upper CI"
	lsmdif    label = "LS Mean Difference"
	lsmdlci   label = "Lower 95% CI of LS Mean Difference"
	lsmduci   label = "Upper 95% CI of LS Mean Difference"
	pval      label = "LS Mean Difference P-Value";
  length analyid analydsc pop studyid subgroup category trt  $ 200;
  format trtn best8.;
  length timeptw covstr $ 200;
  format timeptwn n lsmean lsmeanse lowerci upperci lsmdif lsmdlci lsmduci pval best8.;
  set work.gwms0106_merge
      work.gwcl0403_merge
	  work.gwsp0604_merge
	  work.savant_merge;
  pop = 'EFF';
run;

proc sort data = work.mmrmnrss
           out = output.mmrmnrss (label='MMRM CFB NRS Spasticity Subgroup' drop = categoryn);
  by analyid studyid subgroup categoryn category descending trtn;
run;

%checklog()
