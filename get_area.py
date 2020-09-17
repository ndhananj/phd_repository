
from modes import *

def get_area(xvgfile,area_string=r"$Bottleneck Area (\AA^2)$"):
    A=calc_triangulated_surface_area(xvgfile)
    plot_area(A,area_string=area_string)

if __name__ == '__main__':
    xvgfile = sys.argv[1]
    area_string = sys.argv[2] if len(sys.argv)>2 \
       else r"$Bottleneck Area (\AA^2)$"
    get_area(xvgfile,area_string)
