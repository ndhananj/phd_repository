from pymol import cmd, stored

################################################################################
# pymol visualization functions
################################################################################

@cmd.extend
def show_only_non_waters(water_res="SOL"):
    cmd.hide(representation="everything", selection="all")
    cmd.select(name="non_water", selection="not(resn "+water_res+")")
    cmd.show(representation="sticks", selection="non_water")

@cmd.extend
def load_sparse_traj(traj,obj=None,state=1,interval=400):
    if obj==None:
        [obj,ext]=traj.split('.')
    my_cmd = 'load_traj '+traj+', '+obj+', '+str(state)+\
       ', interval='+str(interval)
    cmd.do(my_cmd)
