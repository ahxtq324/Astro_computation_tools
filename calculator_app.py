import webbrowser
import socket
import subprocess
import os
from flask import Flask, render_template, request
import numpy as np
from scipy.stats import poisson
import scipy.special


def freq(l):
    f = 2.99e10/l
    return f





def vega(F0, mag,dist,efffreq):
    fluxdens = F0*10**(-0.4*mag)*4*np.pi*((dist*3.086e24)**2)*efffreq
    return fluxdens





#app = Flask(__name__)
app = Flask(__name__, template_folder='my_templates')





@app.route('/')
def home():
    return render_template('index.html')





@app.route('/calculate', methods=['GET', 'POST'])
def calculate():
    if request.method == 'POST':
        num1 = float(request.form['num1'])
        num2 = float(request.form['num2'])
        operator = request.form['operator']
        result = None

        if operator == 'add':
            result = num1 + num2
        elif operator == 'subtract':
            result = num1 - num2
        elif operator == 'multiply':
            result = num1 * num2
        elif operator == 'divide':
            result = num1 / num2

        return render_template('calculator.html', result=result)
    return render_template('calculator.html')





@app.route('/opticalmagstolumcalc',methods=['GET', 'POST'])
def opticalmagtolumcalc():
    if request.method == 'POST':
        mag = float(request.form['mag'])
        dist = float(request.form['dist'])
        photometric_system = request.form['Photometric System']
        band = request.form['Band']
        result = None
        
        if photometric_system == "AB":
            if band == "U":
                efffreq = freq(0.36e-4)
            if band == "B":
                efffreq = freq(0.438e-4)
            if band == "V":
                efffreq = freq(0.545e-4)
            if band == "R":
                efffreq = freq(0.641e-4)
            if band == "I":
                efffreq = freq(0.798e-4)
            result = f"{np.power(10, (mag+48.60)/(-2.5))*4*np.pi*((dist*3.086e24)**2)*efffreq:.2e} erg/s "
        
        elif photometric_system == "Vega":
            if band == "U":
                F0 = 1.79e-20 ##values from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
                efffreq = freq(0.36e-4) ## same website as above; in cm
            if band == "B":
                F0 = 4.063e-20 ##values from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
                efffreq = freq(0.438e-4) ## same website as above; in cm
            if band == "V":
                F0 = 3.636e-20 ##values from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
                efffreq = freq(0.545-4) ## same website as above; in cm
            if band == "R":
                F0 = 3.064e-20 ##values from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
                efffreq = freq(0.641e-4) ## same website as above; in cm 
            if band == "I":
                F0 = 2.416e-20 ##values from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
                efffreq = freq(0.798e-4) ## same website as above; in cm
            result = f"{vega(F0,mag,dist,efffreq):0.3e} erg/s"
        return render_template('opticalmagtolumcalc.html', result=result)
    return render_template('opticalmagtolumcalc.html')





@app.route('/fluxtofluxdensity',methods=['GET', 'POST'])
def fluxtofluxdensity():
    if request.method == 'POST':
        fl_ergcm2s = float(request.form['F'])               #r'Flux (F; erg/cm$^2$/s) to be converted to flux density (F$_{\nu}$)'])
        fl_perr = float(request.form['F_perr'])             #r'Positive error on $\delta$F$_+$ (erg/cm$^2$/s)'])
        fl_nerr = float(request.form['F_nerr'])             #r'Negative error on $\delta$F$_-$ (erg/cm$^2$/s)'])
        ene_low = float(request.form['ene_low'])            #Lower energy bound (keV)'])
        ene_high = float(request.form['ene_high'])          #Higher energy bound (keV)'])
        ene_fldens = float(request.form['ene_fldens'])      #'Energy at which F$_{\nu}$ needs to be calculated (keV)']
        gamma = float(request.form['gamma'])                #Best-fitting X-ray photon index']
        result = None

        if gamma != 2:
            conv    = 2.42e17
            nu_low  = ene_low * conv
            nu_high = ene_high * conv
            nu_dens = ene_fldens * conv
            beta    = gamma - 1 #1.1  ##F = nu^(-beta)

            fl_ergcm2shz    = (fl_ergcm2s * (1-beta) * (nu_dens) ** (-beta)) / (nu_high ** (1.-beta) - nu_low ** (1.-beta))
            fl_mjy           = fl_ergcm2shz * 1e26
            fl_mjy_perr       = (fl_perr/fl_ergcm2s) * fl_mjy
            fl_mjy_nerr       = (fl_nerr/fl_ergcm2s) * fl_mjy
            print ('{%.2e}'%fl_mjy + r'$^{+%.2e}$'%(fl_mjy_perr) + r'$_{-%.2e}$'%(fl_mjy_nerr))
            result = f"{fl_mjy:.1e} (+{fl_mjy_perr:.1e},-{fl_mjy_nerr:.1e}) mJy" #'{%.2e}'%fl_mjy + r'$^{+%.2e}$'%(fl_mjy_perr) + r'$_{-%.2e}$'%(fl_mjy_nerr)
        elif gamma == 2:
            print ('The photon index shoud not be equal to 2. If nothing else, try to use 2.0001')
            result = 'The photon index shoud not be equal to 2. If nothing else, try to use 2.0001 or 1.9999'
        return render_template('fluxtofluxdensity.html', result=result)
    return render_template('fluxtofluxdensity.html')





@app.route('/fluxdensitytoflux',methods=['GET', 'POST'])
def fluxdensitytoflux():
    if request.method == 'POST':
        fl_nu_mjy = float(request.form['F_nu'])                  #r'Flux (F; erg/cm$^2$/s) to be converted to flux density (F$_{\nu}$)'])
        flnu_perr = float(request.form['Fnu_perr'])             #r'Positive error on $\delta$F$_+$ (erg/cm$^2$/s)'])
        flnu_nerr = float(request.form['Fnu_nerr'])             #r'Negative error on $\delta$F$_-$ (erg/cm$^2$/s)'])
        ene_low = float(request.form['ene_low'])            #Lower energy bound (keV)'])
        ene_high = float(request.form['ene_high'])          #Higher energy bound (keV)'])
        ene_fldens = float(request.form['ene_fldens'])      #'Energy at which F$_{\nu}$ needs to be calculated (keV)']
        spidx = float(request.form['spidx'])                #Best-fitting X-ray photon index']
        result = None

        if spidx != 1:
            conv    = 2.42e17
            nu_low  = ene_low * conv
            nu_high = ene_high * conv
            nu_dens = ene_fldens * conv
            beta    = spidx #1.1  ##F = nu^(-beta)

            fl_ergcm2sHz          =  fl_nu_mjy/1e26
            fl_ergcm2s_1            = (fl_ergcm2sHz / (nu_dens * conv)**(-1.*beta)) * ((nu_high**(1.-beta) - nu_low**(1.-beta)) / (1.-beta))
            fl_ergcm2s_perr_1       = (flnu_perr/fl_nu_mjy) * fl_ergcm2s_1
            fl_ergcm2s_nerr_1       = (flnu_nerr/fl_nu_mjy) * fl_ergcm2s_1
            print ('{%.2e}'%fl_ergcm2s_1 + r'$^{+%.2e}$'%(fl_ergcm2s_perr_1) + r'$_{-%.2e}$'%(fl_ergcm2s_nerr_1))
            result = f"{fl_ergcm2s_1:.2e} (+{fl_ergcm2s_perr_1:.2e},-{fl_ergcm2s_nerr_1:.2e}) erg/cm2/s" #'{%.2e}'%fl_mjy + r'$^{+%.2e}$'%(fl_mjy_perr) + r'$_{-%.2e}$'%(fl_mjy_nerr)
        elif spidx == 1:
            print ('The spectral index shoud not be equal to 1. If nothing else, try to use 2.0001')
            result = 'The spectral index shoud not be equal to 1. If nothing else, try to use 1.0001 or 0.9999'
        return render_template('fluxdensitytoflux.html', result=result)
    return render_template('fluxdensitytoflux.html')





@app.route('/significanceofdetection',methods=['GET', 'POST'])
def significanceofdetection():
    if request.method == 'POST':
        expected_bckcounts = float(request.form['bkg_counts'])              #]'What is the expected background counts in the source region? ')
        totalcounts        = float(request.form['total_counts'])          #input('What is the total count # in the source area? ')
        result = None

        if totalcounts - expected_bckcounts >= 1:
            poissonprob = 1. - poisson.cdf(int(totalcounts),expected_bckcounts) +  poisson.pmf(int(totalcounts),expected_bckcounts)
            sigma = scipy.special.erfinv(1. - 2.*poissonprob) * np.sqrt(2)
            result = f"Source counts above background is {int(float(totalcounts) - expected_bckcounts)}. \nThe (Gaussian equivalent) significance of detection is {sigma:.2f} σ."
        else:
            print ('The source is not detected')
            result = result = f"Source counts above background is {int(float(totalcounts) - expected_bckcounts)} \n The source is not detected."
        return render_template('significanceofdetection.html', result=result)
    return render_template('significanceofdetection.html')





@app.route('/threesigmaupperlimit',methods=['GET', 'POST'])
def threesigmaupperlimit():
    if request.method == 'POST':
        expected_bckcounts = float(request.form['bkg_counts'])              
        exposure_time      = float(request.form['expo_time'])
        sigma_level      = float(request.form['sigma_level'])  
        result = None

        confidence_level = scipy.special.erf(sigma_level/np.sqrt(2))
        Total_counts_upper = poisson.ppf(confidence_level,expected_bckcounts)

        result = f"Total counts that need to be observed for a 3σ (gaussian equivalent) detection is {Total_counts_upper}. \n3σ count-rate upper-limit is {Total_counts_upper/(float(exposure_time)*1000.):.2e} c/s."
        return render_template('threesigmaupperlimit.html', result=result)
    return render_template('threesigmaupperlimit.html')





if __name__ == '__main__':
    port = 5008

    # Check if the port is available
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('localhost', port))
    if result == 0:
        print(f"Port {port} is already in use. Attempting to terminate the existing process...")
        
        # Kill the existing process on the port
        try:
            os.system(f"fuser -k {port}/tcp")
            #os.system(f"kill -9 $(lsof -t -i:{port})")
            #subprocess.run(['lsof', '-ti', f'tcp:{port}', '|', 'xargs', 'kill'])
            print("Existing process terminated successfully.")
        except subprocess.CalledProcessError:
            print("Failed to terminate the existing process.")
        sock.close()
        
        
    # Specify the path to the Google Chrome executable on Mac
    chrome_path = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'

    # Open the website in Google Chrome
    webbrowser.register('chrome', None, webbrowser.BackgroundBrowser(chrome_path))
    webbrowser.get('chrome').open_new_tab(f"http://localhost:{port}")
    app.run(debug=True,port=port,use_reloader=False)

