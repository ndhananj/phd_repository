################################################################################
# This script is to be called from pymol to load the visualization state
#    showing the comparison of Complex I between 6rfr and 4hea
# Usage (from pymol) : run <this_file_name>
# Example : run load_visualization_state.py
################################################################################

from pymol import cmd, stored

cmd.space("cmyk")
cmd.bg_color(color="white")

cmd.load('4hea_vis_res.pdb')
cmd.load('6rfr_vis_res.pdb')
cmd.load('refined_cav.pdb')
cmd.load('6rfr_UQ9_realigned.pdb')

cmd.hide(representation="everything", selection="all")
cmd.show(representation="mesh", selection="refined_cav")
cmd.show(representation="sticks", selection="4hea_vis_res")
cmd.show(representation="spheres", selection="6rfr_UQ9_realigned")
cmd.show(representation="sticks", selection="6rfr_vis_res")

cmd.do('color cyan, 6rfr_UQ9_realigned and name C*')
cmd.do('color cyan, 6rfr_vis_res and name C*')
cmd.do('color pink, 4hea_vis_res and name c*')
cmd.do('color pink, refined_cav')

cmd.zoom("all", complete=1 )
cmd.rotate("z", angle=90)
cmd.rotate("x", angle=90)

pov_file = open('pymol_scene_no_hydrogens_tyr_marker_set.pov','w')
for part in cmd.get_povray():
    pov_file.write(part)
pov_file.close()

cmd.quit()
