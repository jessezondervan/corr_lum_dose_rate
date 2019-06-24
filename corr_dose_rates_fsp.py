#### feldspar
### calculate the age based on a model for water content in terrace sediments
##need file with input data, see example csv file
fil = open('C:/Users/jzondervan/OneDrive - University of Plymouth/Risoe 03-06-2019/csv_corr_ages_IR50.csv')

#set the following variables: water content error (fractional), grain attenuation of beta dose rates(fractional), 
#internal dose rate (Gy/ky) and error (fractional), internal dose rate of feldspar (Gy/ky) and error (fractional), 
# beta calibration error (fractional)

wc_err = 0.04
grain_att = 0.887355
internal = 0.1
internal_se = 0.05
internal_fsp = 0.857
internal_fsp_err = 0.041
beta_calib = 0.02

####################################################################################################
# Function written by Jesse R. Zondervan - Updated : 24/06/19
####################################################################################################

fil.readline()

sample_numbers = []

age_samples = []
age_ends = []

wc_sats = []
wc_wets = []

gamma_drys = []
g_ses = []
beta_drys = []
b_ses = []
conglom_corrs = []

external_gammas = []
ext_ses = []

Des = []
De_ses = []

for line in fil:
    sample_numbers.append(line.split(',')[0])
    age_samples.append(float(line.split(',')[1]))
    age_ends.append(float(line.split(',')[2]))
    wc_sats.append(float(line.split(',')[3]))
    wc_wets.append(float(line.split(',')[4]))
    gamma_drys.append(float(line.split(',')[5]))
    g_ses.append(float(line.split(',')[6]))
    beta_drys.append(float(line.split(',')[7]))
    b_ses.append(float(line.split(',')[8]))
    conglom_corrs.append(float(line.split(',')[9]))
    external_gammas.append(float(line.split(',')[10]))
    ext_ses.append(float(line.split(',')[11]))
    Des.append(float(line.split(',')[12]))
    De_ses.append(float(line.split(',')[13]))
            


############################### start of code ############
########### Jesse Zondervan 14 June 2019 #################
import numpy as np

for i in range(0,len(age_samples)):
    
    age_sample = age_samples[i]
    age_end = age_ends[i]

    wc_sat = wc_sats[i]/100
    wc_wet = wc_wets[i]/100

    gamma_dry = gamma_drys[i]
    g_se = g_ses[i]
    beta_dry = beta_drys[i]
    b_se = b_ses[i]
    conglom_corr = conglom_corrs[i]
    beta_dry = beta_dry *conglom_corr
    b_se = b_se*conglom_corr

    external_gamma = external_gammas[i]
    ext_se = ext_ses[i]

    De = Des[i]
    De_se = De_ses[i]


    time_sat = age_sample - age_end
    time_wet = age_end

    wc_att_gamma_sat = 1/(1+(1.14*wc_sat))
    wc_att_gamma_wet = 1/(1+(1.14*wc_wet))
    wc_att_beta_sat = 1/(1+(1.25*wc_sat))
    wc_att_beta_wet = 1/(1+(1.25*wc_wet))

    gamma_sat = gamma_dry*wc_att_gamma_sat + external_gamma
    gamma_wet = gamma_dry*wc_att_gamma_wet + external_gamma
    beta_sat = beta_dry*wc_att_beta_sat * grain_att 
    beta_wet = beta_dry*wc_att_beta_wet * grain_att

    tot_rate_sat = gamma_sat + beta_sat + internal + internal_fsp
    tot_rate_wet = gamma_wet + beta_wet + internal + internal_fsp

    weighted_rate = (tot_rate_sat * (time_sat/age_sample)) + (tot_rate_wet * (time_wet/age_sample))

    age_corr = De/weighted_rate

    while abs(age_corr - age_sample) > 2:
        age_sample = age_corr
        time_sat = age_sample - age_end
        time_wet = age_end

        weighted_rate = (tot_rate_sat * (time_sat/age_sample)) + (tot_rate_wet * (time_wet/age_sample))

        age_corr = De/weighted_rate
    
    gamma_sat_se = np.sqrt((g_se*wc_att_gamma_sat)**2+ext_se**2)
    beta_sat_se = b_se*wc_att_beta_sat*grain_att
    gamma_wet_se = np.sqrt((g_se*wc_att_gamma_wet)**2+ext_se**2)
    beta_wet_se = b_se*wc_att_beta_wet*grain_att
    tot_rate_sat_se = np.sqrt(gamma_sat_se**2+beta_sat_se**2+internal_se**2+internal_fsp_err**2)
    tot_rate_wet_se = np.sqrt(gamma_wet_se**2+beta_wet_se**2+internal_se**2+internal_fsp_err**2)


    #now calculate the water content error
    wc_sat = wc_sat+wc_err
    wc_wet = wc_wet+wc_err

    wc_att_gamma_sat = 1/(1+(1.14*wc_sat))
    wc_att_gamma_wet = 1/(1+(1.14*wc_wet))
    wc_att_beta_sat = 1/(1+(1.25*wc_sat))
    wc_att_beta_wet = 1/(1+(1.25*wc_wet))

    gamma_sat = gamma_dry*wc_att_gamma_sat + external_gamma
    gamma_wet = gamma_dry*wc_att_gamma_wet + external_gamma
    beta_sat = beta_dry*wc_att_beta_sat * grain_att 
    beta_wet = beta_dry*wc_att_beta_wet * grain_att

    tot_rate_sat = gamma_sat + beta_sat + internal + internal_fsp
    tot_rate_wet = gamma_wet + beta_wet + internal + internal_fsp

    weighted_rate_xtra = (tot_rate_sat * (time_sat/age_sample)) + (tot_rate_wet * (time_wet/age_sample))

    age_corr_xtra = De/weighted_rate_xtra

    while abs(age_corr_xtra - age_sample) > 5:
        age_sample = age_corr_xtra
        time_sat = age_sample - age_end
        time_wet = age_end

        weighted_rate_xtra = (tot_rate_sat * (time_sat/age_sample)) + (tot_rate_wet * (time_wet/age_sample))

        age_corr_xtra = De/weighted_rate_xtra
    
    wat_cont_err = abs(age_corr_xtra-age_corr)/age_corr

    weighted_rate_se = np.sqrt(tot_rate_sat_se**2+tot_rate_wet_se**2)

    age_corr_se_prec_only = np.sqrt((De_se/De)**2+weighted_rate_se**2)*age_corr
    age_corr_se_all = np.sqrt((De_se/De)**2+weighted_rate_se**2+beta_calib**2+wat_cont_err**2)*age_corr
    
    print (sample_numbers[i])
    print ('the corrected age is %.f +- %.f kyr (prec only is %.f kyr)' % (age_corr, age_corr_se_all, age_corr_se_prec_only))
    print ('the corrected weighted rate is %.2f +- %.2f Gy/kyr' % (weighted_rate, weighted_rate_se))
