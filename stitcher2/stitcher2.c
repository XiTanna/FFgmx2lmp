/* Objective: 
   To build topology of an MD system by combining the topology of each "type" of molecules 
   
   To run:
   ./stitcher data.system Ntype m1.itp Nm1 molecular? m2.itp Nm2 molecular? ... allMolecule.pos
   
   Input:
   1. topolgy of each molecule in the Lammps_tpl format, molA.itp
   2. number of each type of molecules 
   3. whether molecule is described by the "molecular" type in Lammps
   4. the coordinate file of all in unit of angstroms (.pos file); 
      the order of molecules must EXACTLY match the command line input
   
   format of .pos file
   1st line: number of atoms in the file
   after that, each line has three numbers (x, y, z)
   
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

	
   Assumptions:

   1. atom, bond, angle, dihedrals of different molecule types, even same, will be counted as 
      different(e.g., CH2 atoms in octane and hextane are the same, but will be treated as 
      different) to simplify this code

   2. if angle/bond/dihedral/improper is missing, then its heading should NOT be included!!!


*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "types/state.h"
#include "types/smalloc.h"

void mmic(FILE *name);
int ** CreateMatrixInt(int rows,int cols); double ** CreateMatrix(int rows,int cols);
int reportEachMolecule(t_mol *mol,int n_molType);
int totAtomBeforeThisMolecule(t_mol *mol,int i,int k,int *molAmount);

int main(int argc, char *argv[])
{
  int i,j,k,n1,n2,n3,count,nCheck;
  double r1,r2,r3;
  char temp[128], line[1024];
  char delims[] = " ";	
  char *word    = NULL;	
  
  t_mol	*mol;
  FILE *sysFile  = fopen(argv[1],"w");
  FILE *posFile  = fopen(argv[argc-1],"r");
  int n_molTypes = atoi(argv[2]);
  
  int *molAmount;           // number of molecules in the system for each molecule type
  int totAtomType,totBondType,totAngType,totDihType,totImpType;
  int totAtoms,totBonds,totAngs,totDihs,totImps;
  int totMassType;
  int atomID,molID,bondID,angleID,dihID,impID;
  int molecularType;        // whether we will take the entire system as "molecular" type
  int bondTypeBeforeThis,atmTypeBeforeThis,angTypeBeforeThis,dihTypeBeforeThis,impTypeBeforeThis;
  int atmBeforeThisMolecule;
  printf("to run: ./stitcher outFile molTypes m1.itp n_m1 molecular? m2.itp n_m2 molecular? ... all.pos \n");
  printf("\n1. the .itp file for each molecule type must begin with strictly xxx atoms\n");

  printf("There are %3d types of molecules inside your system\n",n_molTypes);
  if((argc-4)%3!=0){
    printf("./stitcher outFile molTypes m1.itp n_m1 molecular? m2.itp n_m2 molecular? ... all.pos \n");
    exit(0);
  }

  // read from separate files
  snew(molAmount,n_molTypes); 
  snew(mol,n_molTypes);      // create pointer aiming at each molecule type
  for(i=0;i<n_molTypes;i++){
    mol[i].fp      = fopen(argv[3+i*3],"r");
    molAmount[i]   = atoi(argv[3+i*3+1]);
    mol[i].molecular=atoi(argv[3+i*3+2]);

    // 1. preambles
    fscanf(mol[i].fp,"%d %s",&mol[i].nAtoms,temp);
    fscanf(mol[i].fp,"%d %s",&mol[i].nBonds,temp);
    fscanf(mol[i].fp,"%d %s",&mol[i].nAngles,temp);
    fscanf(mol[i].fp,"%d %s",&mol[i].nDihedral,temp);
    fscanf(mol[i].fp,"%d %s",&mol[i].nImp,temp);

    fscanf(mol[i].fp,"%d %s %s",&mol[i].atomTypes,temp,line);
    fscanf(mol[i].fp,"%d %s %s",&mol[i].bondTypes,temp,line);
    fscanf(mol[i].fp,"%d %s %s",&mol[i].angleTypes,temp,line);
    fscanf(mol[i].fp,"%d %s %s",&mol[i].dihTypes,temp,line);
    fscanf(mol[i].fp,"%d %s %s",&mol[i].impTypes,temp,line);

    // 2. mass information
    fscanf(mol[i].fp,"%s",temp);
    snew(mol[i].massAtom,mol[i].atomTypes);

    for(j=0;j<mol[i].atomTypes;j++)
      fscanf(mol[i].fp,"%d %lf",&k,&mol[i].massAtom[j]);
    
    // 3. atom informaiton. will need to check molecular format though!
    snew(mol[i].typAtom,mol[i].nAtoms);
    snew(mol[i].chgAtom,mol[i].nAtoms);

    // from first line data, determine whehter this is a "molecular" type
    fscanf(mol[i].fp,"%s",temp);
    for(j=0;j<mol[i].nAtoms;j++){
      if(mol[i].molecular==1){
	fscanf(mol[i].fp,"%d %d %d %lf %lf %lf",&n1,&n2,&mol[i].typAtom[j],&r1,&r2,&r3);
	mol[i].chgAtom[j] = 0.0;}
      else
	fscanf(mol[i].fp,"%d %d %d %lf %lf %lf %lf",&n1,&n2,&mol[i].typAtom[j],&mol[i].chgAtom[j],&r1,&r2,&r3);
    }				
    
    // 4. bonds information
    if(mol[i].nBonds>0){
      fscanf(mol[i].fp,"%s",temp);
      mol[i].bonds = CreateMatrixInt(mol[i].nBonds,3);
      for(j=0;j<mol[i].nBonds;j++)
	fscanf(mol[i].fp,"%d %d %d %d",&n1,&mol[i].bonds[j][0],&mol[i].bonds[j][1],&mol[i].bonds[j][2]);
    }

    // 5. angles information
    if(mol[i].nAngles>0){
      fscanf(mol[i].fp,"%s",temp);
      mol[i].angles = CreateMatrixInt(mol[i].nAngles,4);
      for(j=0;j<mol[i].nAngles;j++)
	fscanf(mol[i].fp,"%d %d %d %d %d",&n1,&mol[i].angles[j][0],&mol[i].angles[j][1],&mol[i].angles[j][2],&mol[i].angles[j][3]);
    }

    // 6. dihedral information
    if(mol[i].nDihedral>0){
      fscanf(mol[i].fp,"%s",temp);
      mol[i].dih = CreateMatrixInt(mol[i].nDihedral,5);
      for(j=0;j<mol[i].nDihedral;j++)
	fscanf(mol[i].fp,"%d %d %d %d %d %d",&n1,&mol[i].dih[j][0],&mol[i].dih[j][1],&mol[i].dih[j][2],&mol[i].dih[j][3],&mol[i].dih[j][4]);
    }

    // 7. improper information
    if(mol[i].nImp>0){
      fscanf(mol[i].fp,"%s",temp);
      mol[i].imp = CreateMatrixInt(mol[i].nImp,5);
      for(j=0;j<mol[i].nImp;j++)
	fscanf(mol[i].fp,"%d %d %d %d %d %d",&n1,&mol[i].imp[j][0],&mol[i].imp[j][1],&mol[i].imp[j][2],&mol[i].imp[j][3],&mol[i].imp[j][4]);
    }
    }

	// position information
	int totAtomCount = 0;
	for(i=0;i<n_molTypes;i++)
		totAtomCount += mol[i].nAtoms*molAmount[i];
	double **atomPosition = CreateMatrix(totAtomCount,3);
	
	fscanf(posFile,"%d",&nCheck);
	if(nCheck != totAtomCount){
		printf("the number of atoms computed from command line is %8d while that in the position input is %8d",totAtomCount,nCheck);
		printf("check your file ...\n");
		exit(0);
	}
	for(j=0;j<nCheck;j++)
		fscanf(posFile,"%lf %lf %lf",&atomPosition[j][0],&atomPosition[j][1],&atomPosition[j][2]);	

  // print onto screen to check whether data are read as expected
  reportEachMolecule(mol,n_molTypes);
  
  // generate topology of system
  // 1. preamble
  fprintf(sysFile,"LAMMPS Description\n\n");
  
  totAtoms = totBonds = totAngs = totDihs  = totImps  = 0;
  for(i=0;i<n_molTypes;i++){
    totAtoms += mol[i].nAtoms * molAmount[i];
    totBonds += mol[i].nBonds * molAmount[i];
    totAngs  += mol[i].nAngles* molAmount[i];
    totDihs  += mol[i].nDihedral*molAmount[i];
    totImps  += mol[i].nImp   * molAmount[i];
  }
  fprintf(sysFile,"%8d  atoms\n",totAtoms);
  fprintf(sysFile,"%8d  bonds\n",totBonds);
  fprintf(sysFile,"%8d  angles\n",totAngs);
  fprintf(sysFile,"%8d  dihedrals\n",totDihs);
  fprintf(sysFile,"%8d  impropers\n\n",totImps);

  totAtomType = totBondType = totAngType = totDihType = totImpType = 0;
  for(i=0;i<n_molTypes;i++){
    totAtomType += mol[i].atomTypes;
    totBondType += mol[i].bondTypes;
    totAngType  += mol[i].angleTypes;
    totDihType  += mol[i].dihTypes;
    totImpType  += mol[i].impTypes;
  }
  fprintf(sysFile,"%8d  atom types\n",totAtomType);
  fprintf(sysFile,"%8d  bond types\n",totBondType);
  fprintf(sysFile,"%8d  angle types\n",totAngType);
  fprintf(sysFile,"%8d  dihedral types\n",totDihType);
  fprintf(sysFile,"%8d  improper types\n\n",totImpType);

  fprintf(sysFile,"num_for_x_low\t number_for_x_high\t xlo xhi\n");
  fprintf(sysFile,"num_for_y_low\t number_for_y_high\t ylo yhi\n");
  fprintf(sysFile,"num_for_z_low\t number_for_z_high\t zlo zhi\n");
  fprintf(sysFile,"\n");

  fprintf(sysFile,"Masses\n\n");
  totMassType = 1;
  for(i=0;i<n_molTypes;i++){
    for(j=0;j<mol[i].atomTypes;j++){
      fprintf(sysFile,"%4d %8.3f \n",totMassType,mol[i].massAtom[j]);
      totMassType++;
    }
  }
  fprintf(sysFile,"\n");

  // 2. atoms
  fprintf(sysFile,"Atoms\n\n");
  atomID = molID = 1;
  molecularType =  1;
  for(i=0;i<n_molTypes;i++) molecularType *= mol[i].molecular; // "molecular" only if all constitutents are "molecular"
  for(i=0;i<n_molTypes;i++){
    atmTypeBeforeThis = 0;
    for(j=0;j<i;j++)
      atmTypeBeforeThis += mol[j].atomTypes;

    for(k=0;k<molAmount[i];k++){
      for(j=0;j<mol[i].nAtoms;j++){
		if(molecularType==1)  // must use atomID-1 becasue c's numbering starts from 0
			fprintf(sysFile,"%9d\t%9d\t%9d\t%12.6f\t%12.6f\t%12.6f\n",atomID,molID,mol[i].typAtom[j]+atmTypeBeforeThis,
                                                           atomPosition[atomID-1][0],atomPosition[atomID-1][1],atomPosition[atomID-1][2]);
		else
			fprintf(sysFile,"%9d\t%9d\t%9d\t%8.6f\t%12.6f\t%12.6f\t%12.6f\n",atomID,molID,mol[i].typAtom[j]+atmTypeBeforeThis,mol[i].chgAtom[j],
                                                           atomPosition[atomID-1][0],atomPosition[atomID-1][1],atomPosition[atomID-1][2]);
		atomID++;
		// printf("atomID = %d\n",atomID);
      }
      molID++;
    }    
  }

  // 3. bonds
  if(totBonds>0)
    fprintf(sysFile,"\nBonds\n\n");
  bondID = 1;
  for(i=0;i<n_molTypes;i++){
    bondTypeBeforeThis = 0;
    for(j=0;j<i;j++) bondTypeBeforeThis += mol[j].bondTypes;
    for(k=0;k<molAmount[i];k++){
      // number of atoms before the k-th instance of i-th type molecules
      atmBeforeThisMolecule = totAtomBeforeThisMolecule(mol,i,k,molAmount);
      for(j=0;j<mol[i].nBonds;j++){
	fprintf(sysFile,"%9d %9d %9d %9d\n",bondID,mol[i].bonds[j][0]+bondTypeBeforeThis,
		mol[i].bonds[j][1]+atmBeforeThisMolecule,mol[i].bonds[j][2]+atmBeforeThisMolecule);
	bondID++;
      }
    }
  }

  // 4. angles
  if(totAngs>0)
    fprintf(sysFile,"\nAngles\n\n");
  angleID = 1;
  for(i=0;i<n_molTypes;i++){
    angTypeBeforeThis = 0;
    for(j=0;j<i;j++) angTypeBeforeThis += mol[j].angleTypes;
    for(k=0;k<molAmount[i];k++){
      // number of atoms before the k-th instance of i-th type molecules
      atmBeforeThisMolecule = totAtomBeforeThisMolecule(mol,i,k,molAmount);
      for(j=0;j<mol[i].nAngles;j++){
	fprintf(sysFile,"%9d %9d %9d %9d %9d\n",angleID,mol[i].angles[j][0]+angTypeBeforeThis,
		mol[i].angles[j][1]+atmBeforeThisMolecule,mol[i].angles[j][2]+atmBeforeThisMolecule,
		mol[i].angles[j][3]+atmBeforeThisMolecule);
	angleID++;
      }
    }
  }

  // 5. dihedrals
  if(totDihs>0)
    fprintf(sysFile,"\nDihedrals\n\n");
  dihID = 1;
  for(i=0;i<n_molTypes;i++){
    dihTypeBeforeThis = 0;
    for(j=0;j<i;j++) dihTypeBeforeThis += mol[j].dihTypes;
    for(k=0;k<molAmount[i];k++){
      // number of atoms before the k-th instance of i-th type molecules
      atmBeforeThisMolecule = totAtomBeforeThisMolecule(mol,i,k,molAmount);
      for(j=0;j<mol[i].nDihedral;j++){
	fprintf(sysFile,"%9d %9d %9d %9d %9d %9d\n",dihID,mol[i].dih[j][0]+dihTypeBeforeThis,
		mol[i].dih[j][1]+atmBeforeThisMolecule,mol[i].dih[j][2]+atmBeforeThisMolecule,
		mol[i].dih[j][3]+atmBeforeThisMolecule,mol[i].dih[j][4]+atmBeforeThisMolecule);
	dihID++;
      }
    }
  }

  // 6. impropers
  if(totImps>0)
    fprintf(sysFile,"\nImpropers\n\n");
  impID = 1;
  for(i=0;i<n_molTypes;i++){
    impTypeBeforeThis = 0;
    for(j=0;j<i;j++) impTypeBeforeThis += mol[j].impTypes;
    for(k=0;k<molAmount[i];k++){
      // number of atoms before the k-th instance of i-th type molecules
      atmBeforeThisMolecule = totAtomBeforeThisMolecule(mol,i,k,molAmount);
      for(j=0;j<mol[i].nImp;j++){
	fprintf(sysFile,"%9d %9d %9d %9d %9d %9d\n",impID,mol[i].imp[j][0]+impTypeBeforeThis,
		mol[i].imp[j][1]+atmBeforeThisMolecule,mol[i].imp[j][2]+atmBeforeThisMolecule,
		mol[i].imp[j][3]+atmBeforeThisMolecule,mol[i].imp[j][4]+atmBeforeThisMolecule);
	impID++;
      }
    }
  }

  fclose(sysFile);

  return 0;
}

int totAtomBeforeThisMolecule(t_mol *mol,int i,int k,int *molAmount)
{
  int count = 0;
  int m = 0;
  for(m=0;m<i;m++)
    count += mol[m].nAtoms * molAmount[m];

  if(k>0)
    count += mol[i].nAtoms * k;

  return count;
}

int reportEachMolecule(t_mol *mol,int n_molType)
{
  int i,j;
  for(i=0;i<n_molType;i++){

    printf("\n===================== Molecule %3d ======================\n",i+1);
    printf("%d atoms\n",mol[i].nAtoms);
    printf("%d bonds\n",mol[i].nBonds);
    printf("%d angles\n",mol[i].nAngles);
    printf("%d dihedrals\n",mol[i].nDihedral);
    printf("%d impropers\n",mol[i].nImp);
    
    printf("%d atom types\n",mol[i].atomTypes);
    printf("%d bond types\n",mol[i].bondTypes);
    printf("%d angle types\n",mol[i].angleTypes);
    printf("%d dihedral types\n",mol[i].dihTypes);
    printf("%d improper types\n",mol[i].impTypes);
    
    printf("Masses \n");
    for(j=0;j<mol[i].atomTypes;j++)
      printf("%d %lf\n",j+1,mol[i].massAtom[j]);
    
    printf("Atoms\n");
    for(j=0;j<mol[i].nAtoms;j++)
      printf("%d %d %d %lf \n",j+1,1,mol[i].typAtom[j],mol[i].chgAtom[j]);
    
    if(mol[i].nBonds>0){
      printf("Bonds\n");
      for(j=0;j<mol[i].nBonds;j++)
	printf("%d %d %d %d\n",j+1,mol[i].bonds[j][0],mol[i].bonds[j][1],mol[i].bonds[j][2]);
    }

    if(mol[i].nAngles>0){
      printf("Angles\n");
      for(j=0;j<mol[i].nAngles;j++)
	printf("%d %d %d %d %d\n",j+1,mol[i].angles[j][0],mol[i].angles[j][1],mol[i].angles[j][2],mol[i].angles[j][3]);
    }

    if(mol[i].nDihedral>0){
      printf("Dihedrals\n");
      for(j=0;j<mol[i].nDihedral;j++)
	printf("%d %d %d %d %d %d\n",j+1,mol[i].dih[j][0],mol[i].dih[j][1],mol[i].dih[j][2],mol[i].dih[j][3],mol[i].dih[j][4]);
    }

    if(mol[i].nImp>0){
      printf("Impropers\n");
      for(j=0;j<mol[i].nImp;j++)
	printf("%d %d %d %d %d %d\n",j+1,mol[i].imp[j][0],mol[i].imp[j][1],mol[i].imp[j][2],mol[i].imp[j][3],mol[i].imp[j][4]);		
    }
  }
  
	return 0;
}
// ----- do not touch these codes below ---------------
// to skip the comment in input files
void mmic(FILE *name)
{
 char d;
 do
 {
 fscanf(name,"%c",&d);
 }
 while(d!='*');
 return;
}
int ** CreateMatrixInt(int rows,int cols)
{
   int  i;
   int  **m;

   m = calloc((unsigned int) rows,sizeof(int *));
   for (i=0; i < rows; i++) {
      m[i] = calloc((unsigned int) cols,sizeof(int));
   }

   return m;
}

double ** CreateMatrix(int rows,int cols)
{
   int  i;
   double  **m;

   m = calloc((unsigned int) rows,sizeof(double *));
   for (i=0; i < rows; i++) {
      m[i] = calloc((unsigned int) cols,sizeof(double));
   }

   return m;
}
