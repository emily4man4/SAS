proc datasets library = work kill nolist;
run;
quit;

data work.data;
  set work.adsk;
run;

proc lifetest data = work.data method = km atrisk plot = (survival);
  time aval*cnsr(1); *in this case, 1 indicates a censored observation;
  *strata sex; *add in this step when you want to present by sex;
  ods output ProductLimitEstimates = work.estimates;
run;

proc lifetest data = work.data method = km atrisk alpha = 0.05 timelist = (0 to 48 by 4); *Update timelist option here;
  time aval*cnsr(1); 
  *strata sex;
  ods output ProductLimitEstimates = work.timelist;
run;

*calculate survival at each time point, for plotting graph;
data work.estimates1 (drop = temp censor survival);
  set work.estimates;
  surv=survival;
  retain temp;
  if surv~=. then temp=surv;
  else surv=temp;
  if censor=1 then cens=aval;
run;

data work.estimates2 (keep = aval cens surv);
  set work.estimates1 (where = (cens ne .));
run;

proc sort data = work.estimates2 out = work.estimates3 nodupkey;
  by aval cens;
run;

* Remove observations that do not need to be in the table;
data work.estimates4;
 set work.estimates1 (where = (numberatrisk ne .));
run;

data work.estimates5 (keep = aval cens surv ord1);
  merge work.estimates3 
        work.estimates4;
  by aval cens;
  ord1=1;
run;

data work.timelist1 (keep = timelist numberatrisk ord1);
  set work.timelist;
  ord1=2;
run;

data work.km;
  format aval timelist 8.;
  attrib _all_ label='';

  set work.estimates5 
      work.timelist1;
run;

title1 "Kaplan-Meier Plot";
proc sgplot data = work.km;
  step x = aval y = surv / name = "s"; *add option, group = sex, here when presenting by sex;
  scatter x = cens y = surv / markerattrs = (symbol = plus)
                              name = "c";
  scatter x = cens y = surv / markerattrs = (symbol = plus); *add option, group = sex, here when presenting by sex;
  xaxistable numberatrisk / x = timelist
						    title = "Number at Risk"
						    location = outside
                            position = bottom; *add option, class = sex, here when presenting by sex;
  xaxis label = "Time (months)";
  yaxis label = "Survival Distribution Function"
        values = (0 to 1 by 0.25);
  keylegend "c" / location = inside 
                  position = topright;
  *keylegend "s"; *add this back in when presenting by sex;
run;
