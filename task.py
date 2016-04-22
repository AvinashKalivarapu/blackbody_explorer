#Author : Avinash Kalivarapu

from bokeh.io import vform,hplot,vplot
from bokeh.plotting import Figure, output_file, show
from bokeh.models import HoverTool,CustomJS, Slider,ColumnDataSource
from bokeh.models.widgets import Panel, Tabs, Dropdown,Select
import numpy as np

h = 6.6244e-34 # Planck constant in Js
c = 2.998e8     # speed of light in m/s
k = 1.38e-23    #Boltzmann constant in J/K
pi = 3.14159     
temp = 5500
jj = 1e-9
a = 1.0
b= 1.0
D = 3.0856e17 #As D = 10pc

output_file("task.html")

gv = []
gv.append(c)
gv.append(h)
gv.append(k)
gv.append(jj)
gv.append(pi)
gv.append(b)
gv.append(D)
for i in range(1,494):
    gv.append(i)
#print gv[0]
#print gv[1]
#print gv[2]
#print gv[3]
#print len(gv)

def bb_flux(wv,t):
    a = 2.0*h*(c)**2/(wv*jj)**5
    b = np.exp(h*c/(k*t*wv*jj))-1.0
    return a/b

ux_pass = np.arange(300,425,5)
uy_pass = [0.00,0.016,0.068,0.167,0.287,0.423,0.560,0.673,0.772,0.841,0.905,0.943,0.981,0.993,1.000,0.989,0.916,0.804,0.625,0.423,0.238,0.114,0.051,0.019,0.000]
uz_pass = bb_flux(ux_pass,temp) * uy_pass

bx_pass = np.arange(360,570,10)
by_pass = [0.0,0.030,0.134,0.567,0.920,0.978,1.000,0.978,0.935,0.853,0.740,0.640,0.536,0.424,0.325,0.235,0.150,0.095,0.043,0.009,0.0]
bz_pass = bb_flux(bx_pass,temp) * by_pass

vx_pass = np.arange(470,710,10)
vy_pass = [0.000,0.030,0.163,0.458,0.780,0.967,1.000,0.973,0.898,0.792,0.684,0.574,0.461,0.359,0.270,0.197,0.135,0.081,0.045,0.025,0.017,0.013,0.009,0.000]
vz_pass = bb_flux(vx_pass,temp) * vy_pass


rx_pass = np.arange(550,910,10)
ry_pass = [0.000,0.23,0.74,0.91,0.98,1.000,0.98,0.96,0.93,0.90,0.86,0.81,0.78,0.72,0.67,0.61,0.56,0.51,0.46,0.40,0.35,0.301,0.277,0.232,0.187,0.14,0.123,0.091,0.079,0.055,0.03,0.025,0.0147,0.0119,0.009,0.00]
rz_pass = bb_flux(rx_pass,temp) * ry_pass
#rz_pass = []

def bb_filter(mean,stde,x):
    d = x-mean
    a = (-(d)**2 ) /2.0*(stde)**2
    b = np.sqrt(2*pi)*(stde*jj)
    c = np.exp(a)/b
    return c

def bb_filter_my(m,s,x):
    return 1/(s * jj * np.sqrt(2*pi)) * np.exp(-((x-m))**2 / (2*s**2))


def Filtters(mean,FWHM):
    start = mean-FWHM
    end = mean+FWHM
    deviation = FWHM/2.355


'''#.................................Filters_My_Computations.......................................................

#For Violet 
#mean = 551nm
#FWHM = 88nm
#standard_dev = 88/2.335

v_mean = 551.0
v_start = 551.0-88.0
v_end = 551.0+88.0
v_dev = 88.0/2.355

v_x = np.linspace(v_start,v_end,1000)
v_y = bb_filter_my(v_mean,v_dev,v_x)

#For Blue.................
#mean = 445nm
#FWHM = 94nm
#standard_dev = 94.0/2.335

b_mean = 445.0
b_start = 445.0-94.0
b_end = 445.0+94.0
b_dev = 94.0/2.355

b_x = np.linspace(b_start,b_end,1000)
b_y = bb_filter_my(b_mean,b_dev,b_x)

#For ultraviolet............
#mean = 365nm
#FWHM = 66nm
#standard_dev = 66.0/2.355
#print u_dev
aub = np.sqrt(2.0*pi) * bb_flux(365,5500)
kkk = aub * jj
#u_dev = 1.0/aub
u_dev = 66.0/2.355
#print u_dev
u_mean = 365.0
u_start = 365.0-66.0
u_end = 365.0+66.0

u_x = np.linspace(u_start,u_end,1000)
#u_y = 1/(u_dev * jj * np.sqrt(2*pi)) * np.exp(-((u_x-u_mean))**2 / (2*u_dev**2))
u_y = bb_filter_my(u_mean,u_dev,u_x)

u_yy = bb_filter(u_mean,u_dev,u_x)
#print u_y
#print u_yy
#print u_y == u_yy
#print bb_flux(365,5500)
testx = 365
testy = 1/(u_dev * np.sqrt(2*pi)) * np.exp(-((testx-u_mean) * jj )**2 / (2*u_dev**2))
#print "................"
#print testy
#print u_x
#print u_y

#For Red..............
#mean = 658nm
#FWHM = 138nm
#standard_dev = 138/2.355

r_mean = 658.0
r_start = 658.0-138.0
r_end = 658.0+138.0
r_dev = 138.0/2.335

r_x = np.linspace(r_start,r_end,1000)
r_y = bb_filter_my(r_mean,r_dev,r_x)
'''

p3 = Figure(plot_width = 500,plot_height = 500)
p3.title = "Plot at 5500K"
p3.yaxis.axis_label = "Flux"
p3.xaxis.axis_label = "Wavelength (nm)"
p3.patch(ux_pass,uz_pass,line_width = 2, line_alpha = 0.1, color = (62,6,148), legend = "Ultraviolet(U)")
p3.patch(bx_pass,bz_pass,line_width = 2, line_alpha = 0.1, color = "blue", legend = "Blue(B)")
p3.patch(vx_pass,vz_pass,line_width = 2, line_alpha = 0.1, color = "violet", legend = "Violet(V)")
p3.patch(rx_pass,rz_pass,line_width = 2, line_alpha = 0.1, color = "red", legend = "Red(R)")



#.....................................FluxVsWavelength Implementation.........................................................

x = np.linspace(50,3000,500)
#print x
y = bb_flux(x,temp)
#print y 
source = ColumnDataSource( data=dict( x=x,y=y,gv=gv))
#source1 = ColumnDataSource( data=dict( x=x,y=y,gv=gv,))
hover = HoverTool( tooltips = [ ("(w,B)", "($x, $y)"), ] )
hover1 = HoverTool( tooltips = [ ("(w,B)", "($x,$y)"), ] )
#p = Figure(plot_width = 500, plot_height = 500, tools = [hover], title="Blackbody" )
p = Figure(plot_width = 500, plot_height = 500, tools = [hover])
#x = np.linspace(50,3000,500)

#y = bb_flux(x,temp)

'''p.line('x','y',source=source, line_width = 2,line_alpha = 0.6)
p.title = "Blackbody curves"
p.yaxis.axis_label = 'Flux (W.sr-1.m-3)'
p.xaxis.axis_label = 'Wavelength (nm)' '''
#p.xaxis.bounds = (400,700)
#tab1 = Panel(child = p, title = "Curves")


for i in range(1,476):

    ux_pass = list(ux_pass)
    uy_pass = list(uy_pass)
    uz_pass = list(uz_pass)
    ux_pass.append(421)
    uy_pass.append(0.00)
    uz_pass.append(0.00)

for i in range(1,480):
    
    bx_pass = list(bx_pass)
    by_pass = list(by_pass)
    bz_pass = list(bz_pass)
    bx_pass.append(561)
    by_pass.append(0.00)
    bz_pass.append(0.00)

for i in range(1,477):
    
    vx_pass = list(vx_pass)
    vy_pass = list(vy_pass)
    vz_pass = list(vz_pass)
    vx_pass.append(701)
    vy_pass.append(0.00)
    vz_pass.append(0.00)

for i in range(1,465):
    
    rx_pass = list(rx_pass)
    ry_pass = list(ry_pass)
    rz_pass = list(rz_pass)
    rx_pass.append(901)
    ry_pass.append(0.00)
    rz_pass.append(0.00)

u21x = [4,2,2,4]
u21y = [0.049,0.049,0,0]
b21x = [7,5,5,7]
b21y = [0.0479,0.0479,0,0]
v21x = [10,8,8,10]
v21y = [0.052,0.052,0,0]
r21x = [13,11,11,13]
r21y = [0.055,0.055,0,0]

for i in range(1,497):
    u21x.append(4)
    u21y.append(0)
    b21x.append(7)
    b21y.append(0)
    v21x.append(10)
    v21y.append(0)
    r21x.append(13)
    r21y.append(0)
    

source1 = ColumnDataSource( data=dict( x=x,y=y,xu=ux_pass,yu=uz_pass,mnu = uy_pass,gv=gv,xb=bx_pass,yb=bz_pass,mnb = by_pass,xv=vx_pass,yv=vz_pass,mnv = vy_pass,xr=rx_pass,yr=rz_pass,mnr = ry_pass,u21x=u21x, u21y=u21y, b21x=b21x, b21y=b21y, v21x=v21x,v21y=v21y,r21x=r21x,r21y=r21y ))

#.................................................JavaScript  Code.......................


callback = CustomJS(args=dict(source=source),code ="""
                var data = source.get('data');
		var f = cb_obj.get('value')
		x = data['x']
		y = data['y']
		gv = data['gv']
		for(i=0; i < x.length; i++)
		{
		y[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(x[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*x[i]*gv[3]))-1.0))
		}
		source.trigger('change');
	""")




callback1 = CustomJS(args=dict(source=source1),code ="""
                var data = source.get('data');
		var f = cb_obj.get('value')
                x = data['x']
                y = data['y']
		xu = data['xu']
		yu = data['yu']
		mnu = data['mnu']
		gv = data['gv']
		xv = data['xv']
		yv = data['yv']
		mnv = data['mnv']
		xb = data['xb']
		yb = data['yb']
		mnb = data['mnb']
		xr = data['xr']
		yr = data['yr']
		mnr = data['mnr']
                u21x = data['u21x']
                u21y = data['u21y']
                b21x = data['b21x']
                b21y = data['b21y']
                v21x = data['v21x']
                v21y = data['v21y']
                r21x = data['r21x']
                r21y = data['r21y']
		for(i=0; i < x.length; i++)
		{
		y[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(x[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*x[i]*gv[3]))-1.0))
		}
		for(i=0; i < xu.length; i++)
		{
		yu[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xu[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xu[i]*gv[3]))-1.0))*mnu[i]
		}
		for(i=0;i < xb.length; i++)
		{
                yb[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xb[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xb[i]*gv[3]))-1.0))*mnv[i]
		}
		for(i=0;i<xv.length; i++)
		{
		yv[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xv[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xv[i]*gv[3]))-1.0))*mnb[i]
		}
		for(i=0;i<xr.length;i++)
		{
		yr[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xr[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xr[i]*gv[3]))-1.0))*mnr[i]
		}
		u21y[0] = 1.0 / ( -2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(370*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*370*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*4.19*1e-20))) 
		u21y[1] = 1.0 / ( -2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(370*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*370*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*4.19*1e-20)))
		b21y[0] = 1.0 / ( -2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(420*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*420*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*6.91*1e-20))) 
		b21y[1] = 1.0 / ( -2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(420*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*420*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*6.91*1e-20)))
		v21y[0] = 1.0 / (-2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(530*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*530*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*3.61*1e-20)))
		v21y[1] = 1.0 / (-2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(530*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*530*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*3.61*1e-20)))
		r21y[0] = 1.0 / (-2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(600*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*600*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*2.19*1e-20)))
		r21y[1] = 1.0 / (-2.5 * Math.log((2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(600*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*600*gv[3]))-1.0))/(4*gv[4]*gv[6]*gv[6]*2.19*1e-20)))
		source.trigger('change');
	""")



#callback1 = CustomJS(args=dict(source=source1),code ="""
#                 var data = source.get('data');
#		 var f = cb_obj.get('value')
#		 x = data['x']
#		 y = data['y']
#		 gv = data['gv']
#       """)
#p.line(u_x,u_y,line_width = 2, line_alpha = 0.6)

#y[i] = (2.0*gv[1]*(gv[0])**2/(x[i]*gv[3])**5)/(Math.exp(gv[1]*gv[0]/(gv[2]*f*x[i]*gv[3]))-1.0)


p.line('x','y',source=source1, line_width = 2,line_alpha = 0.6)
p.title = "Blackbody curves"
p.yaxis.axis_label = 'Flux (W.sr-1.m-3)'
p.xaxis.axis_label = 'Wavelength (nm)'



slider = Slider(start=3000, end = 25000, value = 5500, step = 100, title = "Temperature(in K)", callback= callback1)
#slider1 = Slider(start=3000, end = 25000, value = 5500, step = 100, title = "Temperature(in K)", callback= callback1)


layout = vform(slider,p)
tab1 = Panel(child = layout, title = "Curves" )



#p2 = layout
p2 = Figure(plot_width=580, plot_height=580, tools = [hover1])
p2.line('x','y',source=source1,line_width = 2,line_alpha = 0.6)
#p2.patch(u_x,u_y,line_width = 2, line_alpha = 0.6, color = (62,6,148), legend = "Ultraviolet(U)")
#p2.patch(b_x,b_y,line_width = 2, line_alpha = 0.6, color = "blue", legend = "Blue(B)")
#p2.patch(v_x,v_y,line_width = 2, line_alpha = 0.6, color = "violet", legend = "Violet(V)")
#p2.patch(r_x,r_y,line_width = 2, line_alpha = 0.6, color = "red", legend = "Red(R)")
p2.patch('xu','yu',source=source1,line_width = 2, line_alpha = 0.6, color = (62,6,148), legend = "Ultraviolet(U)")
p2.patch('xb','yb',source=source1,line_width = 2, line_alpha = 0.6, color = "blue", legend = "Blue(B)")
p2.patch('xv','yv',source=source1,line_width = 2, line_alpha = 0.6, color = "violet", legend = "Violet(V)")
p2.patch('xr','yr',source=source1,line_width = 2, line_alpha = 0.6, color = "red", legend = "Red(R)")
p2.title = "Blackbody filters"
p2.yaxis.axis_label = 'Flux (W.sr-1.m-3)'
p2.xaxis.axis_label = 'Wavelength (nm)'

fx = [0]
fy = [0]
p22 = Figure(title="",x_axis_type=None,y_axis_type=None,x_range=(0,20),y_range=(0,20),min_border=0, outline_line_color="black",toolbar_location=None)
p22.circle(fx,fy,color="white",radius = 0.0)
p22.xgrid.grid_line_color = None
p22.ygrid.grid_line_color = None



#U_mag = (bb_flux(370, 5500))/(4*pi*D*D*1.81*1e-23)
U_mag = (bb_flux(370, 5500))/(4*pi*D*D*4.19*1e-20)
u_m = -2.5 * np.log(U_mag)
#print U_mag
#print 1/u_m

#B_mag = (bb_flux(420, 5500))/(4*pi*D*D*4.26*1e-23)
B_mag = (bb_flux(420, 5500))/(4*pi*D*D*6.91*1e-20)
b_m = -2.5 * np.log(B_mag)
#print B_mag
#print 1/b_m

#V_mag = (bb_flux(530, 5500))/(4*pi*D*D*3.64*1e-23)
V_mag = (bb_flux(530, 5500))/(4*pi*D*D*3.61*1e-20)
v_m = -2.5 * np.log(V_mag)
#print V_mag
#print 1/v_m

#R_mag = (bb_flux(600, 5500))/(4*pi*D*D*3.08*1e-23)
R_mag = (bb_flux(600, 5500))/(4*pi*D*D*2.19*1e-20)
r_m = -2.5 * np.log(R_mag)
#print R_mag
#print 1/r_m

#9670
'''
u21x = [2,4,4,2]
u21y = [4,4,0,0]
b21x = [5,7,7,5]
b21y = [6,6,0,0]
v21x = [8,10,10,8]
v21y = [8,8,0,0]
r21x = [11,13,13,11]
r21y = [10,10,0,0]
'''
source21 = ColumnDataSource( data=dict(u21x=u21x, u21y=u21y, b21x=b21x, b21y=b21y, v21x=v21x,v21y=v21y,r21x=r21x,r21y=r21y,gv=gv ))

'''callback21 = CustomJS(args=dict(source=source1),code ="""
                var data = source.get('data');
		var f = cb_obj.get('value')
                x = data['x']
                y = data['y']
		xu = data['xu']
		yu = data['yu']
		mnu = data['mnu']
		gv = data['gv']
		xv = data['xv']
		yv = data['yv']
		mnv = data['mnv']
		xb = data['xb']
		yb = data['yb']
		mnb = data['mnb']
		xr = data['xr']
		yr = data['yr']
		mnr = data['mnr']
		for(i=0; i < x.length; i++)
		{
		y[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(x[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*x[i]*gv[3]))-1.0))
		}
		for(i=0; i < xu.length; i++)
		{
		yu[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xu[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xu[i]*gv[3]))-1.0))*mnu[i]
		}
		for(i=0;i < xb.length; i++)
		{
                yb[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xb[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xb[i]*gv[3]))-1.0))*mnv[i]
		}
		for(i=0;i<xv.length; i++)
		{
		yv[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xv[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xv[i]*gv[3]))-1.0))*mnb[i]
		}
		for(i=0;i<xr.length;i++)
		{
		yr[i] = (2.0*gv[1]*Math.pow(gv[0],2)/(Math.pow(xr[i]*gv[3],5))/(Math.exp((gv[1]*gv[0])/(gv[2]*f*xr[i]*gv[3]))-1.0))*mnr[i]
		}
		source.trigger('change');
	""")'''




p21 = Figure(title="relative strengths at temperatures",x_axis_type=None,y_axis_type=None,x_range=(0,20),y_range=(0,0.5),min_border=0, outline_line_color="black",toolbar_location=None)
p21.patch('u21x','u21y',source=source1,line_width = 2,line_alpha = 0.6,color = (62,6,148),legend = "U")
p21.patch('b21x','b21y',source=source1,line_width = 2,line_alpha = 0.6,color = "blue",legend = "B")
p21.patch('v21x','v21y',source=source1,line_width = 2,line_alpha = 0.6,color = "violet",legend = "V")
p21.patch('r21x','r21y',source=source1,line_width = 2,line_alpha = 0.6,color = "red",legend = "R")

p21.xgrid.grid_line_color = None
p21.ygrid.grid_line_color = None


menu = [("U", "item_1"), ("V", "item_2"), ("B","item_3"), ("R", "item_4")]
dropdown = Dropdown(label="U", menu=menu)
select = Select( value="U - U", options=["U - U", "U - B", "U - V", "U - R", "B - U", "B - B", "B - V", "B - R", "V - U", "V - B", "V - V","V - R", "R - U", "R - B", "R - V", "R - R"])

layout_p2 = vform(slider,p2)
layout_p22 = vform(select,p21)
layout1 = hplot(layout_p2,layout_p22)

tab2 = Panel(child=layout1, title="filters")

tab3 = Panel(child=p3, title = "Detailed Plot")


#.........................................Transmission passbands........................
#ux_pass = np.arange(300,425,5)
#uy_pass = [0.00,0.016,0.068,0.167,0.287,0.423,0.560,0.673,0.772,0.841,0.905,0.943,0.981,0.993,1.000,0.989,0.916,0.804,0.625,0.423,0.238,0.114,0.051,0.019,0.000]

#bx_pass = np.arange(360,570,10)
#by_pass = [0.0,0.030,0.134,0.567,0.920,0.978,1.000,0.978,0.935,0.853,0.740,0.640,0.536,0.424,0.325,0.235,0.150,0.095,0.043,0.009,0.0]

#vx_pass = np.arange(470,710,10)
#vy_pass = [0.000,0.030,0.163,0.458,0.780,0.967,1.000,0.973,0.898,0.792,0.684,0.574,0.461,0.359,0.270,0.197,0.135,0.081,0.045,0.025,0.017,0.013,0.009,0.000]

#rx_pass = [550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,800,850,900] 
#ry_pass = [0.000,0.23,0.74,0.91,0.98,1.000,0.98,0.96,0.93,0.90,0.86,0.81,0.78,0.72,0.67,0.61,0.56,0.51,0.46,0.40,0.35,0.14,0.03,0.00]

source_u = ColumnDataSource( data=dict( x=ux_pass,y=uy_pass))
source_b = ColumnDataSource( data=dict( x=bx_pass,y=by_pass))
source_v = ColumnDataSource( data=dict( x=vx_pass,y=vy_pass))
source_r = ColumnDataSource( data=dict( x=rx_pass,y=ry_pass))
source_pass = ColumnDataSource( data=dict( x=[ux_pass,bx_pass,vx_pass,rx_pass], y=[uy_pass,by_pass,vy_pass,ry_pass], desc = ['U','B','V','R']))


hover_u = HoverTool( tooltips = [ ("(wavelen,Trans)", "($x, $y)"), ] )
hover_b = HoverTool( tooltips = [ ("(wavelen,Trans)", "($x, $y)"),("Description","Red passband"), ])
#hover_pass = HoverTool( tooltips = [ ("(wavelength,Trans)","(@$x,@$y)"),("Description","@desc"))


p4 = Figure(plot_width=500, plot_height=500, tools = [hover_u])
p4.line(ux_pass,uy_pass,line_width = 2,line_alpha = 0.6, source = source_u, color = (62, 6, 148), legend = "Ultraviolet(U)")
p4.line(bx_pass,by_pass,line_width = 2,line_alpha = 0.6, source = source_b, color = "blue", legend = "Blue(B)")
p4.line(vx_pass,vy_pass,line_width = 2,line_alpha = 0.6, source = source_v, color = "violet", legend = "Violet(V)")
p4.line(rx_pass,ry_pass,line_width = 2,line_alpha = 0.6, source = source_r, color = "red", legend = "Red(R)")

p4.title = "Besell approximations UBVR passbands"
p4.yaxis.axis_label = 'Transmission'
p4.xaxis.axis_label = 'Wavelength (nm)'

tab4 = Panel(child=p4, title = "Transmissions")


tabs = Tabs(tabs=[ tab1, tab2, tab3, tab4 ])

#show(tabs)


#show(layout)
show(tabs)
#show(p)

