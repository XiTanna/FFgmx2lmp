#include <stdlib.h>
#include <stdio.h>

typedef struct
{
  FILE    *fp;          /* file containing topology in Lammps format   */
  FILE    *fp2;         /* file containing the coodinates of each atom */

  int     molecular;    /* ==1 means molecular, i.e., no charge        */
  int	  atomTypes;    /* total type of atoms in ...                  */
  double  *massAtom;    /* mass of each type of atoms                  */

  int     nAtoms;       /* number of atoms in this molecule            */
  double  *chgAtom;     /* charge of each atom in this molecule        */
  int     *typAtom;     /* type of atom -- for each atom in molecule   */
  
  int     nBonds;       /* number of bonds in ...                      */
  int 	  bondTypes;    /* total type of bonds in ...                  */
  int     **bonds;      /* bond[i][0]: type, bond[i][1-2]: atom_1-2    */
  
  int 	  nAngles;      /* number of angles in ...                     */
  int 	  angleTypes;   /* total type of angles in ...                 */
  int     **angles;      /* angle[i][0]: type, angle[i][1-3]: atom_1-3  */
  
  int 	  nDihedral;    /* number of dihedrals in ...                  */
  int 	  dihTypes;     /* total type of dihedrals in ...              */
  int     **dih;        /* dih[i][0]: type, dih[i][1-4]: atom_1-4      */
    
  int 	  nImp;         /* number of impropers in ...                  */
  int 	  impTypes;     /* total type of impropers in ...              */
  int     **imp;        /* imp[i][0]:type, imp[1-4]: atom_1-4          */
  
  double **pos;         /* coordinate of each atom in the system       */

} t_mol;

/* pos file of coordinates   
   line 1: total atom number
   line 2-N: coordinate of each atom in angstrom unit
*/

/*
topology in Lammps_tpl format:
   
    33016 atoms                       [atoms in a single molecule]
    59192 bonds                       [bonds in a single ..      ]
    43560 angles
    65524 dihedrals
        0 impropers

        3 atom types
        3 bond types
        2 angle types
        3 dihedral types
        0 improper types

	Masses
	1 xxx
	2 xxx
	...

	Atoms
	
	1	1	1	1077.989	4.970	3.400
	2	1	2	1077.702	5.538	4.787
	...

	Bonds
	1	xx	yy yy
	...

	Angles
	1	xx yy yy yy 
	...
	
	Dihedrals
	1	xx yy yy yy yy
	...
	
	Impropers
	1	xx yy yy yy yy
*/
