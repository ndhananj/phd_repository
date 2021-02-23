
from pymol import cmd, stored

cmd.hide(representation="everything", selection="all")
cmd.show(representation="sticks", selection="all")

cmd.do('spectrum b, white_red')

#cmd.zoom("all", complete=1 )
#cmd.rotate("z", angle=90)
#cmd.rotate("x", angle=90)
