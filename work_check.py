
from gmx_file_processing import *

def checkTimeMatch(dist_df,force_df):
    return np.all(0==(dist_df[['Time (ps)']].to_numpy()\
       -force_df[['Time (ps)']].to_numpy()))

def calcWork1D_Trap(x,f):
    dx=(x[1:]-x[:-1]).flatten()
    f_avg = 0.5*(f[1:]+f[:-1]).flatten()
    x_avg = 0.5*(x[1:]+x[:-1]).flatten()
    integrand = -(f_avg*dx)
    data_dict = {'x':x_avg, 'dW':integrand, 'work':np.cumsum(integrand)}
    return pd.DataFrame(data=data_dict)


if __name__ == '__main__':
    dist_file = sys.argv[1]
    force_file = sys.argv[2]
    dist_df = read_xvg(dist_file)['data']
    force_df = read_xvg(force_file)['data']
    if(checkTimeMatch(dist_df,force_df)):
       print("The Time axes match")
       x = dist_df[['Position (nm)']].to_numpy()
       f = force_df[['Force (kJ/mol/nm)']].to_numpy()
       work_df = calcWork1D_Trap(x, f)
       fig = plt.figure()
       ax1 = fig.add_subplot(211)
       ax1.scatter('x','dW', data=work_df)
       plt.ylabel('dW (kJ/mol)')
       ax2 = fig.add_subplot(212)
       ax2.scatter('x','work', data=work_df)
       plt.xlabel('x (nm)')
       plt.ylabel('work (kJ/mol)')
       plt.show()
    else:
       print("The Time axes don't match")
