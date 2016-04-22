#include "sg/superglue.hpp"

double ** Amatrix;
double ** Bmatrix;
double ** Cmatrix;

int dim=9000;
int bl_dim=100;

int nu_blck=dim/bl_dim;

void print();

template<typename Options>

struct MyHandle : public HandleBase<Options>{
	size_t i, j;
	void set(size_t i_, size_t j_){
		i=i_;
		j=j_;
	}
	size_t geti(){return i;}
	size_t getj(){return j;}
};

struct Options : public DefaultOptions<Options>{
	typedef MyHandle<Options> HandleType;
};


void blck_Identification(int a1, int a2, int b1, int b2){

	int row_start=a1/bl_dim;
	int row_end=a2/bl_dim-1;
	if (a2%bl_dim) ++row_end;

	int col_start=b1/bl_dim;
	int col_end=b2/bl_dim-1;
	if (b2%bl_dim) ++col_end;
}

struct Mul_Task : public Task<Options> {

	int a1, a2, b1, b2, row_start, row_end, col_start, col_end,
	    a1p, a2p, b1p, b2p, row_start_p, row_end_p, col_start_p, col_end_p;
	Mul_Task(int a1_, int a2_, int b1_, int b2_, int row_start_, int row_end_, int col_start_, int col_end_,
			int a1p_, int a2p_, int b1p_, int b2p_, int row_start_p_, int row_end_p_, int col_start_p_, int col_end_p_,
			Handle<Options> **HA, Handle<Options> **HB, Handle<Options> **HC)
		: a1(a1_), a2(a2_), b1(b1_), b2(b2_), row_start(row_start_), row_end(row_end_), col_start(col_start_), col_end(col_end_),
		a1p(a1p_), a2p(a2p_), b1p(b1p_), b2p(b2p_), row_start_p(row_start_p_), row_end_p(row_end_p_), col_start_p(col_start_p_), col_end_p(col_end_p_)

	{
		for (int i=row_start; i<=row_end; ++i){
			for(int j=col_start; j<=col_end; ++j){
				register_access(ReadWriteAdd::write, HA[i][j]);
				register_access(ReadWriteAdd::write, HB[i][j]);
				register_access(ReadWriteAdd::write, HC[i][j]);

			}
		}
	}

	void run(){
		for (int i= a1; i<a2; ++i){
			for (int j=b1; j<b2; ++j){
				for (int k=a1p; k<a2p; ++k){
					Cmatrix[i][j]=Cmatrix[i][j]+Amatrix[i][k]*Bmatrix[k][j];
				}
			}
		}


		for (int i= a1; i<a2; ++i){
			for (int j=b1; j<b2; ++j){
				Amatrix[i][j]=Cmatrix[i][j];
				Cmatrix[i][j]=0;
			}
		}

	}

};

struct Mul_Task_Col : public Task<Options> {

	int a1, a2, b1, b2, row_start, row_end, col_start, col_end,
	    a1p, a2p, b1p, b2p, shift_b1, shift_b2, row_start_p, row_end_p, col_start_p, col_end_p, j_index;
	Mul_Task_Col(int a1_, int a2_, int b1_, int b2_, int row_start_, int row_end_, int col_start_, int col_end_,
			int a1p_, int a2p_, int b1p_, int b2p_, int shift_b1_, int shift_b2_, int row_start_p_, int row_end_p_, int col_start_p_, int col_end_p_, int j_index_,
			Handle<Options> **HA, Handle<Options> **HB, Handle<Options> **HC)
		: a1(a1_), a2(a2_), b1(b1_), b2(b2_), row_start(row_start_), row_end(row_end_), col_start(col_start_), col_end(col_end_),
		a1p(a1p_), a2p(a2p_), b1p(b1p_), b2p(b2p_), shift_b1(shift_b1_), shift_b2(shift_b2_), row_start_p(row_start_p_), row_end_p(row_end_p_), col_start_p(col_start_p_), col_end_p(col_end_p_), j_index(j_index_)

	{
		for (int i=row_start; i<=row_end; ++i){
			for(int j=col_start; j<=col_end; ++j){
				register_access(ReadWriteAdd::write, HA[i][j]);
				register_access(ReadWriteAdd::write, HB[i][j]);
				register_access(ReadWriteAdd::write, HC[i][j]);

			}
		}
	}

	void run(){



		int jp=b1+(j_index*bl_dim);



		std::cout<<""<<"a1  "<<a1<<'\t'<<"a2  "<<a2<<'\t'<<"b1  "<<shift_b1<<'\t'<<"b2  "<<shift_b2<<
			std::endl<<"a1p "<<a1p<<'\t'<<"a2p "<<a2p<<'\t'<<"b1p "<<b1p<<'\t'<<"b2p "<<b2p<<std::endl;


		for (int i= a1; i<a2; ++i){

			for (int j=b1p, jp=b1+(j_index*bl_dim); j<b2p; ++j, ++jp){

				for (int k=shift_b1; k<shift_b2; ++k){

					int kp=a1p;


					Cmatrix[i][jp]=Cmatrix[i][jp]+Amatrix[i][k]*Bmatrix[kp][j];

					kp++;
				}
			}
		}
	}
};


void row_mul(int a1, int a2, int b1, int b2,
		int a1p, int a2p, int b1p, int b2p, Handle <Options> **HA, Handle <Options> **HB, Handle <Options> **HC){

	int row_start=a1/bl_dim;
	int row_end=a2/bl_dim-1;
	if (a2%bl_dim) ++row_end;

	int col_start=b1/bl_dim;
	int col_end=b2/bl_dim-1;
	if (b2%bl_dim) ++col_end;

	int row_start_p=a1p/bl_dim;
	int row_end_p=a2p/bl_dim-1;
	if (a2p%bl_dim) ++row_end_p;

	int col_start_p=b1p/bl_dim;
	int col_end_p=b2p/bl_dim-1;
	if (b2p%bl_dim) ++col_end_p;


	SuperGlue<Options> sg;

	int step=(a1/bl_dim)+1;
	int na1=a1;
	int na2=step*bl_dim;

	// Row broadcaster
	for (int i=1 ; na1<a2; ){
		sg.submit(new Mul_Task(na1, na2, b1, b2, row_start, row_end, col_start, col_end,
					a1p, a2p, b1p, b2p, row_start_p, row_end_p, col_start_p, col_end_p,
					HA, HB, HC));


		na1=(i+1)*bl_dim;
		na2=(i+2)*bl_dim;
		if (na2>a2) na2=a2;
		i++;
	}

	// Wait for all tasks to finish
	sg.barrier();



}

void col_mul(int a1, int a2, int b1, int b2,
		int a1p, int a2p, int b1p, int b2p, Handle <Options> **HA, Handle <Options> **HB, Handle <Options> **HC){

	int row_start=a1/bl_dim;
	int row_end=a2/bl_dim-1;
	if (a2%bl_dim) ++row_end;

	int col_start=b1/bl_dim;
	int col_end=b2/bl_dim-1;
	if (b2%bl_dim) ++col_end;

	int row_start_p=a1p/bl_dim;
	int row_end_p=a2p/bl_dim-1;
	if (a2p%bl_dim) ++row_end_p;

	int col_start_p=b1p/bl_dim;
	int col_end_p=b2p/bl_dim-1;
	if (b2p%bl_dim) ++col_end_p;

	SuperGlue<Options> sg(1);

	int iter=(b2-b1)/bl_dim;
	int nb1=b1;

	int nb2=b1+bl_dim;
	if (nb2>b2) nb2=b2;

	int na1p=a1p;
	int na2p=a1p+bl_dim;
	if (na2p>a2p) na2p=a2p;

	/// new added
	int nb1p=b1p;
	int nb2p=b1p+bl_dim;
	if (nb2p>b2p) nb2p=b2p;

	int shift_b1=nb1;
	int shift_b2=shift_b1+bl_dim;
	if (shift_b2>b2) shift_b2=b2;


	int sent_task=1;


	/////////////////////////////////

	// finding i,j,k index limit

	int i_lim= row_end - row_start;
	int j_lim= row_end_p - row_start_p;
	int k_lim= col_end_p - col_start_p;

	if (i_lim==0) i_lim++;
	if (j_lim==0) j_lim++;
	if (k_lim==0) k_lim++;

	for ( int i=0; i<i_lim; ++i){
		for (int j=0; j<j_lim; ++j ){
			for (int k=0; k<k_lim; ++k){

				sg.submit(new Mul_Task_Col(a1, a2, nb1, nb2, row_start, row_end, col_start, col_end,
							na1p, na2p, nb1p, nb2p, shift_b1, shift_b2, row_start_p, row_end_p, col_start_p, col_end_p, j, 
							HA, HB, HC));
//				std::cout<<"\n a1 \t"<<a1<<"\n a2 \t"<< a2<<"\n nb1 \t"<< nb1<<"\n nb2 \t"<< nb2<<"\n row_start \t"
//					<< row_start<<"\n row_end \t"<< row_end<<"\n col_start \t"<< col_start<<"\n col_end \t"<< col_end<<
//					"\n na1p \t"<<na1p<<"\n na2p \t"<< na2p<<"\n nb1p \t"<<nb1p<<"\n nb2p \t"<< nb2p<<"\n row_start_p \t"<<
//					row_start_p<<"\n row_end_p \t"<< row_end_p<<"\n col_start_p \t"<< col_start_p<<"\n col_end_p \t"<< col_end_p<<std::endl;


				shift_b1=shift_b2;
				shift_b2+=bl_dim;
				if (shift_b2>b2) shift_b2=b2;

				na1p=na2p;
				na2p+=bl_dim;
				if (na2p>a2p) na2p=a2p;



			}	



			nb1p=nb2p;
			nb2p+=bl_dim;
			if (nb2p>b2p) nb2p=b2p;

			shift_b1=b1;
			shift_b2=b1+bl_dim;
			if (shift_b2>b2) shift_b2=b2;


			nb1=b1;
			nb2=b1+bl_dim;
			if (nb2>b2) nb2=b2;


			na1p=a1p;
			na2p=a1p+bl_dim;
			if (na2p>a2p) na2p=a2p;



		}

		//2std::cout<<"i_lim "<<i_lim<<" "<<"j_lim "<<j_lim<<" "<<"k_lim "<<k_lim<<std::endl;
	}

	// Wait for all tasks to finish
	sg.barrier();

}




void print(){
	for (int i=0; i<dim; ++i){
		for (int j=0; j<dim; ++j){
			std::cout<<Cmatrix[i][j]<<' ';
		}
		std::cout<<std::endl;
	}
}

void set_data(){
	for (int i=0; i<dim; ++i){
		for (int j=0; j<dim; ++j){
			Amatrix[i][j]=0;
		}
	}
}


int main (){
	//std::cout<<"Hello world!"<<std::endl;

	//Matrix Creation
	Amatrix= new double *[(dim)];
	if (Amatrix==NULL){
		std::cout<<"Error in memory allocation";
	}

	for (int i=0; i<dim;++i){
		Amatrix[i]= new double [(dim)];
	}

	//Set Matrix Data
	for (int i=0; i<dim; ++i){
		for (int j=0; j<dim; ++j){
			Amatrix[i][j]=2;
		}
	}

	//Define Handles for the blocks
	Handle <Options> **HA = new Handle <Options> *[nu_blck];
	for (int i=0; i<nu_blck; ++i){
		HA[i]= new Handle <Options> [nu_blck];
	}

	//Matrix Creation
	Bmatrix= new double *[(dim)];
	if (Bmatrix==NULL){
		std::cout<<"Error in memory allocation";
	}

	for (int i=0; i<dim;++i){
		Bmatrix[i]= new double [(dim)];
	}

	//Set 	Matrix Data
	for (int i=0; i<dim; ++i){
		for (int j=0; j<dim; ++j){
			Bmatrix[i][j]=3;
		}
	}

	//Define Handles for the blocks
	Handle <Options> **HB = new Handle <Options> *[nu_blck];
	for (int i=0; i<nu_blck; ++i){
		HB[i]= new Handle <Options> [nu_blck];
	}

	//Matrix Creation
	Cmatrix= new double *[(dim)];
	if (Cmatrix==NULL){
		std::cout<<"Error in memory allocation";
	}

	for (int i=0; i<dim;++i){
		Cmatrix[i]= new double [(dim)];
	}

	//Set Matrix Data
	for (int i=0; i<dim; ++i){
		for (int j=0; j<dim; ++j){
			Cmatrix[i][j]=0;
		}
	}

	//Define Handles for the blocks
	Handle <Options> **HC = new Handle <Options> *[nu_blck];
	for (int i=0; i<nu_blck; ++i){
		HC[i]= new Handle <Options> [nu_blck];
	}
	//*/

	//	blck_Identification(7, 19, 4, 13);

	//	row_mul(7, 107, 4, 13, 4, 13, 4, 13, HA, HB, HC);

	//	col_mul(7, 13, 5, 20, 5, 20, 5, 20, HA, HB, HC);

	//@	col_mul(7, 13, 5, 105, 5, 7995, 5, 7995, HA, HB, HC);

	//#	col_mul(7, 13, 5, 105, 5, 105, 5, 105, HA, HB, HC);

	col_mul(7, 13, 5, 8002, 5, 8002, 5, 8002, HA, HB, HC);

	print();

	delete [] Amatrix; 
	delete [] Bmatrix;
	delete [] Cmatrix;

	delete [] HA;
	delete [] HB;
	delete [] HC;

	return 0;
}

