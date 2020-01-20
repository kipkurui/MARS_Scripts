import numpy
from scipy.optimize import leastsq
import pylab
 
def GetData():
	print("Non-Linear Fitting for Dose-Response Data\n")
	input_x = raw_input("Type X values, comma seperated, then hit return\n")
	input_y = raw_input("Type Y values, comma seperated, then hit return\n")
 
	x = numpy.array(input_x.split(','),dtype=float)
	y = numpy.array(input_y.split(','),dtype=float)
 
	FittingMath(x,y)
	return
 
def FittingMath(x,y):
	equation = lambda variables,x:	(variables[0]+(variables[1]-variables[0])/(1+10**((variables[2]-x))))
	error = lambda variables,x,y: (equation(variables,x)-y)
	variable_guesses = [numpy.min(y),numpy.max(y),numpy.mean(x)]
 
	output = leastsq(error,variable_guesses,args=(x,y),full_output=1 )
	variables = output[0]
	fitinfo = output[2]
	r_squared = residuals(y,fitinfo)
 
	print("The fitted variables are: "+str(variables) )
 
	x_range = numpy.arange(numpy.min(x),numpy.max(x),abs(numpy.max(x)/100))
 
	Grapher(x,y,equation(variables,x_range),x_range,variables[2],r_squared)
	return
 
def residuals(y,fitinfo):
	fit_error = 0
	fit_variance = 0
 
	for i in range(len(fitinfo['fvec'])):
		fit_error += (fitinfo['fvec'][i])**2
		fit_variance += (y[i]-numpy.mean(y))**2
	r_squared = 1 - (fit_error/fit_variance)
	return r_squared
 
def Grapher(x,y,equation,x_range,kd_value,r_squared):
	fig = pylab.figure(1)
	axis = fig.add_subplot(111)
	axis.plot(x,y,'ro',label='data')
	axis.plot(x_range,equation,label='fit')
	axis.set_xlabel('Concentration')
	axis.set_ylabel('Intensity')
	axis.set_title("Dose Response Fit")
	fig.subplots_adjust(bottom=0.15)
	fig.subplots_adjust(left=0.15)
	text = axis.text(0.05, 0.95, 'Kd: %s'% kd_value,transform=axis.transAxes, va='top')
	rtext = axis.text(0.05, 0.90, 'R^2: %s'% r_squared,transform=axis.transAxes, va='top')
	pylab.show()
 
if __name__ == '__main__':
	GetData()