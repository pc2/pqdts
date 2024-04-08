/*
 #  File            : condat_simplexproj.c 
 #
 #  Version History : 1.0, Aug. 15, 2014 
 #
 #  Author          : Laurent Condat, PhD, CNRS research fellow in France.
 #
 #  Description     : This file contains an implementation in the C language
 #                    of algorithms described in the research paper:
 #	
 #                    L. Condat, "Fast Projection onto the Simplex and the
 #                    l1 Ball", preprint Hal-01056171, 2014.
 #
 #                    This implementation comes with no warranty: due to the
 #                    limited number of tests performed, there may remain
 #                    bugs. In case the functions would not do what they are
 #                    supposed to do, please email the author (contact info
 #                    to be found on the web).
 #
 #                    If you use this code or parts of it for any purpose,
 #                    the author asks you to cite the paper above or, in 
 #                    that event, its published version. Please email him if 
 #                    the proposed algorithms were useful for one of your 
 #                    projects, or for any comment or suggestion.
 #
 #  Usage rights    : Copyright Laurent Condat.
 #                    This file is distributed under the terms of the CeCILL
 #                    licence (compatible with the GNU GPL), which can be
 #                    found at the URL "http://www.cecill.info".
 #
 #  This software is governed by the CeCILL license under French law and
 #  abiding by the rules of distribution of free software. You can  use,
 #  modify and or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL :
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
*/

#include <cstdlib>
extern "C"{
void simplexproj_condat_(double* y, double* x, int* plength) {
        int length=plength[0];
	double*	aux = x; //(x==y ? (double*)malloc(length*sizeof(double)) : x);
	double*  aux0=aux;
	int		auxlength=1; 
	int		auxlengthold=-1;	
	double	tau=(*aux=*y)-1.0;
	int 	i=1;
	for (; i<length; i++) 
		if (y[i]>tau) {
			if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))
			<=y[i]-1.0) {
				tau=y[i]-1.0;
				auxlengthold=auxlength-1;
			}
			auxlength++;
		} 
	if (auxlengthold>=0) {
		auxlength-=++auxlengthold;
		aux+=auxlengthold;
		while (--auxlengthold>=0) 
			if (aux0[auxlengthold]>tau) 
				tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
	}
	do {
		auxlengthold=auxlength-1;
		for (i=auxlength=0; i<=auxlengthold; i++)
			if (aux[i]>tau) 
				aux[auxlength++]=aux[i];	
			else 
				tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
	} while (auxlength<=auxlengthold);
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
//	if (x==y) free(aux0);
}
} 
