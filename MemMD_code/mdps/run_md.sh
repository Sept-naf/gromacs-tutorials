gmx editconf -f start_complex.pdb -o box.gro
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
python water_deletor.py
gmx grompp -f ions.mdp -c solv_delW.gro -p topol_delW.top -o ions.tpr
echo -e "SOL\n" | gmx genion -s ions.tpr -o ions.gro -p topol_delW.top -pname NA -nname CL -conc 0.15 -neutral
gmx grompp -f em.mdp -c ions.gro -p topol_delW.top -o em.tpr
gmx mdrun -v -deffnm em

gmx grompp -f nvt.mdp -c em.gro -p topol_delW.top -o nvt.tpr -r em.gro
gmx mdrun -v -deffnm nvt -nt 16

gmx grompp -f npt.mdp -c nvt.gro -p topol_delW.top -o npt.tpr -r nvt.gro
gmx mdrun -v -deffnm npt -nt 16

gmx grompp -f md.mdp -c npt.gro -p topol_delW.top -o md.tpr
gmx mdrun -v -deffnm md -nt 16
