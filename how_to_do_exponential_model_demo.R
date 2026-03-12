
# Example of how to fit an exponential learning model 
# to the median of the response distribution,
# placing boundaries on the rate parameter [the time constant]
#
# (took like 30 minutes on Aaron's laptop; the median regression is relatively slow)

library(brms)

# simulate data with skewed performance distribution:
dat <- do.call('rbind',replicate(10,{
  pars <- c(pStart = runif(1,4,7)
            ,pRate = runif(1,10,50)
            ,pAsym = runif(1,1,3)
  )
  data.frame(
    subID = sample(letters,10,replace=T)
    ,trialNum = 1:200
    ,condition = sample(c('treatment','control'),1)
    ,performance = rnorm(200
                         ,mean = pars['pAsym'] + ( (pars['pStart'])-(pars['pAsym']) )*2^( -(0:199) / (pars['pRate']) ) 
    )^2 # squaring all of the values; SD before squaring was 1
  )
},simplify = F)
)

bForm <- brmsformula(
  performance ~ pAsym + ( pStart-pAsym )*2^( (1 - trialNum) / pRate ) # standard exponential function. the base-2 log could be replaced with base-e or base-10 or whatever, if you want. The `1` in `1-trialNum` should be equal to the smallest value in your time variable
  , pStartLog ~ condition + (1 | subID) # Wilkenson formula for transformed parameter
  , pRateLogis ~ condition + (1 | subID) # Wilkenson formula for transformed parameter
  , pAsymLog ~ condition + (1 | subID) # Wilkenson formula for transformed parameter
  ,nl = T) +
  nlf(pStart ~ pStartLog^2) + # base-2 log transform (ensures positivity; could also do a scaled logistic if you wanted)
  nlf(pRate ~ 5 + (100 - 5) * inv_logit(pRateLogis) ) + # scaled logistic transform (imposes boundaries of 5 and 100)
  nlf(pAsym ~ pAsymLog^2) # base-2 log transform (ensures positivity; could also do a scaled logistic if you wanted)

myModel <- brm(
  bForm
  ,dat
  ,family = asym_laplace(link_quantile = 'identity')
  ,prior = 
    prior(constant(.5),class='quantile') + # ensures estimating median
    prior(normal(0,3),class = 'b', nlpar = 'pStartLog') + # places prior density mostly between 2^-6 and 2^6 [fairly wide]; try `hist(rnorm(1E4,X,X)^2)` to visualize
    prior(normal(0,1.5),class = 'b', nlpar = 'pRateLogis') + # places prior density more-or-less uniformly between our bounds; try `hist(X + (X-X) * plogis(rnorm(1E4,X,X)) )` to visualize alongside the boundaries
    prior(normal(0,3),class = 'b', nlpar = 'pAsymLog') # places prior density mostly between 2^-6 and 2^6 [fairly wide] ; try `hist(rnorm(1E4,X,X)^2)` to visualize
  , chains = 3, cores = 3
  , control = list(adapt_delta = .95)
)  

by_participant_params <- 
  aggregate(cbind(trialNum,performance) ~ subID + condition
            ,myModel$data
            ,median ) ; for(curParam in c('pStart','pRate','pAsym')){
              
              dTmp <- fitted(myModel
                             ,nlpar = curParam
                             ,newdata = by_participant_params)
              by_participant_params[,curParam] <- dTmp[,'Estimate']
              rm(dTmp)
            }

