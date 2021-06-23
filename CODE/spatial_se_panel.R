Spatial_HAC_Panel=function(equation,
                     study_file,
                     range_search,
                     Smoothness=0.5,
                     residual_upper=1000,residual_lower= -1000,
                     opt_method="REML"){
  
  #This computes HAC standard errors for panels. 
  #The datafile MUST have columns called id and year to identify the entries, 
  #and longitude and latitude must be called X and Y.
  #The datafile can be unbalanced by cannot contain missing entries for the regression variables.
  #It will not work if your data are severely unbalanced so that some years have no observations that overlap.
  #The file returns a dataframe called HAC with the regression results, 
  #one called Spatial_Parameters that lists range, structure etc, another called 
  #Annual_Spatial_Params which gives the spatial parameters computed for each year.
  #If the ranges reported in this matrix are above or below the range search values that you specified you should
  #change your range search appropriately.
  #The spatial parameters used are the medians of these values.
  #residual_upper gives a value for truncating large residuals which can distort spatial parameter estimation. 
  
    
study_file=study_file %>% arrange(id,year)
n_years=length(unique(study_file$year))
n_id=length(unique(study_file$id))
set_years=sort(unique(study_file$year))
set_id=sort(unique(study_file$id))

full_study_file=data.frame(id=rep(set_id,each=n_years),year=rep(set_years,n_id))
full_study_file$obs_index=1:nrow(full_study_file)

Coords=study_file %>% group_by(id) %>% summarize(X=min(X,na.rm=T),Y=min(Y,na.rm=T)) %>% select(-id) %>% as.matrix()


ols=lm(as.formula(equation),data=study_file)
Residuals=ols$residuals

##Truncate high values which will mess up spatial parameter ests
Residuals=ifelse(Residuals>residual_upper,residual_upper,Residuals)
Residuals=ifelse(Residuals<residual_lower,residual_lower,Residuals)

study_file$residuals=Residuals
full_study_file=left_join(full_study_file,study_file,by = c("id", "year")) %>% 
  select(id,year,obs_index,residuals)

kappa=Smoothness

test_params=list()
for (i in 1:n_years){
  study_file_yr=study_file %>% filter(year==set_years[i]) %>%  select(X,Y,residuals)%>% na.omit() #%>% 
  locs=study_file_yr %>% select(X,Y) %>% as.matrix()
  
  
  hold_search=data.frame(range_search,lambda=NA,loglik=NA)
  
  for (j in 1:length(range_search)){
    fit_search = Krig(  x = locs,
                        Y = scale(as.vector(study_file_yr$residuals)),
                        Distance = "rdist.earth",
                        Covariance = "Matern", 
                        smoothness = kappa,
                        theta = range_search[j],
                        give.warnings = F
    )
    
    hold_search[j,2:3]=fit_search$lambda.est[6,c(1,5)]
  }
  
  cov_par=hold_search %>% arrange(loglik) %>% slice(1)%>% as.numeric()
  
  study_file_spatial=SpatialPoints(coords=locs)
  proj4string(study_file_spatial)=CRS("+proj=longlat +datum=WGS84")
  nearest=knn2nb(knearneigh(study_file_spatial,k=5,longlat = T))   #k nearest
  nearest=nb2listw(nearest,style="W")
  moran_z=as.vector(moran.test(study_file_yr$residuals,listw=nearest)$statistic)
  cov_par=c(cov_par,moran_z)
  test_params[[i]]=(cov_par)
}

spat_params=ldply(test_params)
names(spat_params)=c("range","lambda","loglik","Moran_z")



##########calc covar mx. Choose median of spatial params for all years
cov_par=as.vector(apply(spat_params,2,median)) %>% as.numeric()

if(cov_par[1]==max(range_search)) warning("Your range search values are too low: try higher ones.")
if(cov_par[1]==min(range_search)) warning("Your range search values are too high: try lower ones.")

Effective_Range_km=1.6*sqrt(8*kappa)*cov_par[1]  #in kilometres
Range=as.numeric(cov_par[1])
Structure=as.numeric(1/(1+cov_par[2]))  #this is weight rho between systematic correlation and spatial noise
loglik=-as.numeric(cov_par[3])
Moran_z=as.numeric(cov_par[4])

spatial_parameters=data.frame(Range,Effective_Range_km,Structure,Moran_z,Smoothness=kappa,
                              loglik)


Cor_Spatial=Structure*fields::Matern(rdist.earth(x1=Coords),
                                  range=Range,smoothness=kappa)+
  diag(nrow(Coords))*(1-Structure)

########calculate temporal autocorrelation matrix Alpha
Alpha=cor(matrix(full_study_file$residuals,byrow=T,ncol=n_years),use="pairwise.complete")

if(is.na(min(Alpha))) stop("Your data are severely unbalanced: some years do not have overlapping observations.")

##Kron product gives spatiotemporal corr
KL=Cor_Spatial %x% Alpha
Index=na.omit(full_study_file)$obs_index   #observations present
KL=KL[Index,Index]

##Calculate HAC standard errors
X=qr.X(ols$qr)
N=nrow(X)
k=ncol(X)
U=ols$res         
V=t(X)%*%(U*KL*U)%*%X/N
sv_X=svd(X)   #invert X'X using SVD
v=sv_X$v
d_2=diag(sv_X$d^(-2))
xx.inv = N *(v %*% d_2 %*% t(v))
Cov=xx.inv%*%V%*%xx.inv/N
hac.se=sqrt(diag(Cov))
hac.t=summary(ols)$coef[,1]/hac.se
hac.p=2*(1-pnorm(abs(hac.t)))
hac=cbind.data.frame(coef=ols$coef,hac.se,hac.t,hac.p)


output=list(HAC=hac,Spatial_Parameters=spatial_parameters,OLS=ols,Residuals=ols$residuals,
            Annual_Spatial_Params=spat_params)
return(output)
}
