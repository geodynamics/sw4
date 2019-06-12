void compute_f( EW& simulation, int nspar, int nmpars, double* xs, int nmpard, double* xm,
		vector<vector<Source*> >& GlobalSources,
		vector<vector<TimeSeries*> >& GlobalTimeSeries,
		vector<vector<TimeSeries*> >& GlobalObservations,
		double& mf, Mopt *mopt );

void compute_f_and_df( EW& simulation, int nspar, int nmpars, double* xs, int nmpard, double* xm,
		       vector<vector<Source*> >& GlobalSources,
		       vector<vector<TimeSeries*> >& GlobalTimeSeries,
		       vector<vector<TimeSeries*> >& GlobalObservations, 
		       double& f, double* dfs, double* dfm, int myrank,
                       Mopt *mopt, int it=-1 );

void linesearch( EW& simulation, vector<vector<Source*> >& GlobalSources,
		 vector<vector<TimeSeries*> >& GlobalTimeSeries, 
		 vector<vector<TimeSeries*> >& GlobalObservations,
		 int nspar, int nmpars, double* xs, int nm, double* xm, double f, double* dfs,
		 double* dfm, double* ps, double* pm, double cgstep, double maxstep, double steptol,
		 double* xsnew, double* xmnew, double& fnew, double* sfs, double* sfm,
		 int myRank, int& retcode, int& nstep_reductions, bool testing, double* dfsnew,
		 double* dfmnew, Mopt* mopt );

void lbfgs( EW& simulation, int nspar, int nmpars, double* xs, 
	    int nmpard, double* xm, 
	    vector<vector<Source*> >& GlobalSources,
	    vector<vector<TimeSeries*> >& GlobalTimeSeries,
	    vector<vector<TimeSeries*> >& GlobalObservations,
	    int myRank, Mopt* mopt );

void nlcg( EW& simulation, int nspar, int nmpars, double* xs,
	   int nmpard, double* xm, 
	   vector<vector<Source*> >& GlobalSources,
	   vector<vector<TimeSeries*> >& GlobalTimeSeries,
	   vector<vector<TimeSeries*> >& GlobalObservations,
	   int myRank, Mopt* mopt );


