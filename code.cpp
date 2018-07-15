#include <mpi.h>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <mpi.h>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <stdio.h>
#include <sys/stat.h>
#include <float.h>
#include <omp.h>
#include <cstdlib>

using namespace std;

typedef struct node
{
	int r;
	int c;
	int valint;
	float valfloat;
}node;


bool structrowcomparator(node a, node b)
{
	// check equal to case
	if(a.r!=b.r)
		return a.r<b.r;
	else 
		return a.c<b.c;
}

bool structcolcomparator(node a, node b)
{
	// check equal to case
	if(a.c!=b.c)
		return a.c<b.c;
	else 
		return a.r<b.r;
}

int intcomparator(const void* a, const void* b)
{
	// check equal to case
	int temp= *(((int*)a) + 1) ;
	int temp1=*(((int*)b) + 1 ) ;
	// cout<<temp<<"ddddddddd"<<endl;

	// cout<<temp[0]<<"    qwertyuio     "<<temp[1] <<endl;

	if(temp<temp1) return -1 ;
	if(temp>temp1) return 1 ;
	return -1 ;

}

int floatcomparator(const void* a, const void* b)
{
	// check equal to case
	float temp= *(((float*)a) + 1) ;
	float temp1=*(((float*)b) + 1) ;
	// cout<<temp<<"ddddddddd"<<endl;

	// cout<<temp[0]<<"    qwertyuio     "<<temp[1] <<endl;

	if(temp<temp1) return -1 ;
	if(temp>temp1) return 1 ;
	return -1 ;

}

// int colcomparator(const void* a, const void* b)
// {
// 	// check equal to case
// 	float* temp=(float*)a;
// 	float* temp1=(float*)b;
// 	// cout<<temp<<"ddddddddd"<<endl;

// 	// cout<<temp[0]<<"    qwertyuio     "<<temp[1] <<endl;

// 	return temp[1]-temp1[1];
// 	// cout<<temp[0]<<endl;
// 	// cout<<"qwertyui"<<endl;


// 	return 1;
// 	// if(a.r!=b.r)
// 	// 	return a.r<b.r;
// 	// else 
// 	// 	return a.c<b.c;
// }


int main(int argc, char *argv[]) 
{
	int num_procs, pid, i=0; 
	float* recv_data;
	int* recv_rc;
	int* recv_cr;
	int* rc_data;
	int* cr_data;
	long fil_sz=0;
	int flag=0;
	// vector< vector<int> > col_inrows;
	// vector< vector<int> > row_incols;

	int* size_ofrows;
	int* size_ofcols;
	float* all_data;
	float* colrecv_data;
	float* colall_data;
	vector<node> all_nodes;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	// MPI_Get_processor_name(proc_name, &name_len);

	float* sorted_cols;
	float* sorted_rows;
	// printf("Hello World from processor %d out of %d, executing on %s\n",pid, num_procs, proc_name);

	// First Process - Input & Output / Distribute Data & Receive Results / Compute
	int num_entries=0;
	int* send_count=(int*)malloc(num_procs*sizeof(int));
	int * displ_arr=(int*)malloc(num_procs*sizeof(int));
	int* colsend_count=(int*)malloc(num_procs*sizeof(int));
	int * coldispl_arr=(int*)malloc(num_procs*sizeof(int));
	

	int rows=0, columns=0;

	if ( pid == 0 ) 
	{
		size_t vecLen;
		char *vec=NULL; char *colI=NULL,*rowI=NULL;

		FILE *inputMat;
		inputMat = fopen(argv[1],"rb");

		// first make a struct and sort in row major form

		fseek (inputMat , 0 , SEEK_END);
		fil_sz = ftell (inputMat);
		rewind (inputMat);
		num_entries=fil_sz/16;
		

		for (i = 0; i < num_entries; i++) 
		{
			// getMatline(&vec,&rowI,&colI,&vecLen,inputMat);

			int* row_n=(int*) malloc(4);
			int* col_n=(int*) malloc(4);
			int* valint=(int*) malloc(4);
			float* valfloat=(float*) malloc(4);

			fread (row_n,4,1,inputMat);
			fread (col_n,4,1,inputMat);
			fread (valint,4,1,inputMat);
			fread (valfloat,4,1,inputMat);

			// assumed in row major form
			node temp;

			temp.valint = (valint[0]) ;
			temp.valfloat = (valfloat[0]);
			temp.c = col_n[0];
			temp.r = row_n[0];

			if(temp.c>columns)
				columns=temp.c;
			if(temp.r>rows)
				rows=temp.r;

			// if (row_n[0] > prevRowI)
			// 	rowptr[j++] = i;  
			// prevRowI = row_n[0];
			all_nodes.push_back(temp);

			free(row_n);
			free(col_n);
			free(valint);
			free(valfloat);
		}
		fclose(inputMat);


		rows++;
		columns++;

		rows=(ceil((rows*1.0f)/num_procs))*num_procs;
		columns=(ceil((columns*1.0f)/num_procs))*num_procs;


		// cout<<rows<<endl<<columns<<endl<<endl<<endl<<endl;

		sort(all_nodes.begin(),all_nodes.end(), structrowcomparator);
		// for(int i=0;i<num_entries;i++)
		// {
		// 	cout<<all_nodes[i].r<<"  "<<all_nodes[i].c<<"  "<<all_nodes[i].valint<<"  "<<all_nodes[i].valfloat<<endl;
		// }

		
		// col_inrows.resize(rows);
		// row_incols.resize(columns);

		size_ofrows=(int*)calloc(rows,sizeof(int));
		size_ofcols=(int*)calloc(columns,sizeof(int));
		
		// for(int i=0;i<num_entries;i++)
		// {
		// 	col_inrows[all_nodes[i].r].push_back(all_nodes[i].c);
		// 	// cout<<all_nodes[i].r<<" "<<all_nodes[i].c<<endl;
		// 	row_incols[all_nodes[i].c].push_back(all_nodes[i].r);
		// 	// cout<<all_nodes[i].c<<" "<<all_nodes[i].r<<endl;
		// }
		// exit(1);

		


		for(int i=0;i<num_entries;i++)
		{
			size_ofrows[all_nodes[i].r]++;
			size_ofcols[all_nodes[i].c]++;
			// cout<<col_inrows[i].size()<<endl;

		}
		// exit(1);

		// cout<<endl<<endl;
		

		// for(int i=0;i<columns;i++)
		// {
			
		// 	// cout<<row_incols[i].size()<<endl;
		// }

		int row_tosend=ceil((rows*1.0f)/num_procs);
		int col_tosend=ceil((columns*1.0f)/num_procs);

		// cout<<"before looper"<<endl;

		// int row_tosend=ceil((rows*1.0f)/num_procs);

		// if(pid==0)
		// {
		all_data=(float*)malloc(2*sizeof(float)*num_entries);
		rc_data=(int*)malloc(2*sizeof(int)*num_entries);
		// }

		for(int i=0;i<num_procs;i++)
		{
			send_count[i]=0;
		}

		for(int i=0;i<rows;i++)
		{
			send_count[i/row_tosend]+=size_ofrows[i]*2;
		}

		displ_arr[0]=0;

		for(int i=1;i<num_procs;i++)
		{
			displ_arr[i]=displ_arr[i-1]+send_count[i-1];
		}


		for(int i=0;i<num_procs;i++)
		{
			colsend_count[i]=0;
		}

		// cout<<"coltosend= "<<col_tosend<<endl<<endl;
		for(int i=0;i<columns;i++)
		{
			colsend_count[i/col_tosend]+=size_ofcols[i]*2;
		}

		coldispl_arr[0]=0;

		for(int i=1;i<num_procs;i++)
		{
			coldispl_arr[i]=coldispl_arr[i-1]+colsend_count[i-1];
		}


	}

		
	MPI_Bcast(send_count,num_procs,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(colsend_count,num_procs,MPI_INT,0,MPI_COMM_WORLD);

	for(int looper=0;looper<4;looper++)
	{
		if(pid==0)
		{
			
			for(int i=0;i<2*num_entries;i+=2)
			{
				int ii = all_nodes[i/2].valint ;
				int* ptr = (int*) (all_data+i) ;
				*ptr = ii ;

				all_data[i+1]=all_nodes[i/2].valfloat;
			}


			for(int i=0;i<2*num_entries;i+=2)
			{
				rc_data[i]=all_nodes[i/2].r;
				rc_data[i+1]=all_nodes[i/2].c;
			}

			
			recv_data=(float*)malloc(send_count[pid]*sizeof(float));
			// recv_data=(float*)malloc(10000000*sizeof(float));
			recv_rc=(int*)malloc(send_count[pid]*sizeof(int));
			// recv_rc=(int*)malloc(10000000*sizeof(int));
			// cout<<"pid 0 "<<endl;
			
			// cout<<"after bcast for 0 "<<pid<<endl;
			MPI_Scatterv(all_data,send_count,displ_arr,MPI_FLOAT,recv_data,send_count[pid],MPI_FLOAT,0,MPI_COMM_WORLD);
			// cout<<"after scatter1 for 0 "<<pid<<endl;
			MPI_Scatterv(rc_data,send_count,displ_arr,MPI_INT,recv_rc,send_count[pid],MPI_INT,0,MPI_COMM_WORLD);
			// cout<<"after scatterv for 0 "<<pid<<endl;
			// for(int i=0;i<send_count[pid];i++)
			// {
			// 	cout<<recv_data[i]<<"   ";
			// }
			// cout<<endl;


		}

		// cout<<"reached after pid=0"<<pid<<endl;


		if(pid!=0)
		{
			// cout<<"pid!=0  "<<pid<<endl;
			// MPI_Bcast(send_count,num_procs,MPI_INT,0,MPI_COMM_WORLD);
			// cout<<"after bcast "<<pid<<endl;
			recv_data=(float*)malloc(send_count[pid]*sizeof(float));
			// recv_data=(float*)malloc(10000000*sizeof(float));
			// float* temp=(float*)malloc(1);
			// int* tempi=(int*)malloc(1);
			MPI_Scatterv(NULL,send_count,NULL,MPI_FLOAT,recv_data,send_count[pid],MPI_FLOAT,0,MPI_COMM_WORLD);
			// cout<<"after scatterv1 "<<pid<<endl;
			recv_rc=(int*)malloc(send_count[pid]*sizeof(int));
			// recv_rc=(int*)malloc(10000000*sizeof(int));
			MPI_Scatterv(NULL,send_count,NULL,MPI_INT,recv_rc,send_count[pid],MPI_INT,0,MPI_COMM_WORLD);
			// cout<<"after scatterv2 "<<pid<<endl;

		}
		// cout<<"reached after pid!=0"<<pid<<endl;

		// MPI_Barrier(MPI_COMM_WORLD);

		// cout<<"reached after 1st barrier"<<endl;



		vector<int> elem_count;
		int rownum=recv_rc[0];
		int row_count=0;
		// calculates num of elements
		// cout<<rownum<<endl;
		for(int i=0;i<send_count[pid];i+=2)
		{
			// cout<<recv_rc[i]<<"ee"<<endl;
			if(recv_rc[i]!=rownum)
			{
				elem_count.push_back(row_count);
				// cout<<row_count<<"yuy"<<endl;
				rownum=recv_rc[i];
				row_count=1;
				if(i==send_count[pid]-1 || i==send_count[pid]-2)
				{
					elem_count.push_back(row_count);
				}

			}
			else
			{
				row_count++;
				if(i==send_count[pid]-1 || i==send_count[pid]-2)
				{
					elem_count.push_back(row_count);
					// cout<<row_count<<"yuy"<<endl;
				}	
			}
		}


		int num_rows=elem_count.size();
		// cout<<num_rows<<endl<<endl<<endl<<endl;
		int* elem_count_prefix=(int*)malloc(sizeof(int)*num_rows);
		elem_count_prefix[0]=0;
		for(int i=1;i<num_rows;i++)
		{
			elem_count_prefix[i]=elem_count_prefix[i-1]+elem_count[i-1]*2;
		}



		#pragma omp parallel for
			for(int i=0;i<num_rows;i++)
			{
				// cout<<"temmmmmmm"<<endl;
				// cout<<elem_count_prefix[i]<<" "<<pid<< endl;
				// cout<<elem_count[i]<<" "<<pid<<endl;
				qsort (recv_data + elem_count_prefix[i], elem_count[i], 2*sizeof(int), floatcomparator);
			}

		int num_objects=0;
		for(int i=0;i<num_procs;i++)
		{
			num_objects+=send_count[i];
		}

		if(pid==0)
		{
			sorted_rows=(float*)malloc(sizeof(float)*num_objects);
		}
		// sorted_rows=(float*)malloc(sizeof(float)*num_objects);

		MPI_Gatherv(recv_data, send_count[pid], MPI_FLOAT,
			sorted_rows, send_count, displ_arr,
	                MPI_FLOAT, 0, MPI_COMM_WORLD);

		// cout<<"reached after 1st gatherv"<<endl;

		// MPI_Barrier(MPI_COMM_WORLD);
		// cout<<"reached after 2nd barrrier"<<endl;

		// printing....................
		// if(pid==0)
		// {
		// 	for(int i=0;i<num_objects;i+=2)
		// 	{
			
		// 		cout<<pid<<": "<<rc_data[i]<<"  "<<rc_data[i+1]<<"->  "<<sorted_rows[i]<<"  "<<sorted_rows[i+1]<< endl;
		// 	}	
		// }

		// here all the data  is in the sorted by rows format in node 0.... now do for columns
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------
		// if(flag==1)
		// 	exit(1);
		// else flag=1;


		if(pid==0)
		{
			
			for(int i=0;i<2*num_entries;i+=2)
			{
				int hhval= *((int*)(sorted_rows+i));
				all_nodes[i/2].valint=hhval;
				all_nodes[i/2].valfloat=sorted_rows[i+1];
			}
			sort(all_nodes.begin(),all_nodes.end(), structcolcomparator);
			

			colrecv_data=(float*)malloc(colsend_count[pid]*sizeof(float));
			colall_data=(float*)malloc(2*sizeof(float)*num_entries);

			for(int i=0;i<2*num_entries;i+=2)
			{
				int ii = all_nodes[i/2].valint ;
				int* ptr = (int*) (colall_data+i+1) ;
				*ptr = ii ;

				colall_data[i]=all_nodes[i/2].valfloat;
				// colall_data[i+1]=all_nodes[i/2].valint;
			}

			cr_data=(int*)malloc(2*sizeof(int)*num_entries);

			for(int i=0;i<2*num_entries;i+=2)
			{
				cr_data[i]=all_nodes[i/2].c;
				cr_data[i+1]=all_nodes[i/2].r;
			}
			
			recv_cr=(int*)malloc(colsend_count[pid]*sizeof(int));
			
			MPI_Scatterv(cr_data,colsend_count,coldispl_arr,MPI_INT,recv_cr,colsend_count[pid],MPI_INT,0,MPI_COMM_WORLD);
			MPI_Scatterv(colall_data,colsend_count,coldispl_arr,MPI_FLOAT,colrecv_data,colsend_count[pid],MPI_FLOAT,0,MPI_COMM_WORLD);

		}

		// cout<<"reached after pid=0 column"<<endl;

		
		if(pid!=0)
		{
			// MPI_Bcast(colsend_count,num_procs,MPI_INT,0,MPI_COMM_WORLD);
			colrecv_data=(float*)malloc(colsend_count[pid]*sizeof(float));


			// float* temp33=(float*)malloc(1);
			// int* tempi33x`=(int*)malloc(1);
			recv_cr=(int*)malloc(colsend_count[pid]*sizeof(int));
			MPI_Scatterv(NULL,colsend_count,NULL,MPI_INT, recv_cr,colsend_count[pid],MPI_INT,0,MPI_COMM_WORLD);
			MPI_Scatterv(NULL,colsend_count,NULL,MPI_FLOAT,colrecv_data,colsend_count[pid],MPI_FLOAT,0,MPI_COMM_WORLD);

		}

		// cout<<"reached after pid!=0 column"<<endl;

		// MPI_Barrier(MPI_COMM_WORLD);

		// cout<<"reached after 3rd barrier"<<endl;

		vector<int> colelem_count;
		int colnum=recv_cr[0];
		int col_count=0;
		// calculates num of elements
		// cout<<rownum<<endl;
		for(int i=0;i<colsend_count[pid];i+=2)
		{
			// cout<<recv_rc[i]<<"ee"<<endl;
			if(recv_cr[i]!=colnum)
			{
				colelem_count.push_back(col_count);
				// cout<<row_count<<"yuy"<<endl;
				colnum=recv_cr[i];
				col_count=1;
				if(i==colsend_count[pid]-1 || i==colsend_count[pid]-2)
				{
					colelem_count.push_back(col_count);
					// cout<<row_count<<endl;
					// cout<<row_count<<"yuy"<<endl;
				}

			}
			else
			{
				col_count++;
				if(i==colsend_count[pid]-1 || i==colsend_count[pid]-2)
				{
					colelem_count.push_back(col_count);
					// cout<<row_count<<"yuy"<<endl;
				}	
			}
		}



		int num_cols=colelem_count.size();
		// cout<<num_rows<<endl<<endl<<endl<<endl;
		int* colelem_count_prefix=(int*)malloc(sizeof(int)*num_cols);
		colelem_count_prefix[0]=0;
		for(int i=1;i<num_cols;i++)
		{
			colelem_count_prefix[i]=colelem_count_prefix[i-1]+colelem_count[i-1]*2;
		}

		#pragma omp parallel for	
		for(int i=0;i<num_cols;i++)
		{
			qsort (colrecv_data + colelem_count_prefix[i], colelem_count[i], 2*sizeof(int), intcomparator);
		}
		

		

		

		int colnum_objects=0;
		for(int i=0;i<num_procs;i++)
		{
			colnum_objects+=colsend_count[i];
		}

		if(pid==0)
		{
			sorted_cols=(float*)malloc(sizeof(float)*colnum_objects);
		}
		

		MPI_Gatherv(colrecv_data, colsend_count[pid], MPI_FLOAT,
			sorted_cols, colsend_count, coldispl_arr,
	                MPI_FLOAT, 0, MPI_COMM_WORLD);

		// cout<<"reached after 2nd gatherv columns"<<endl;


		// MPI_Barrier(MPI_COMM_WORLD);
		// cout<<"reached after 4th barrier"<<endl;

		// if(pid==0)
		// {
		// 	for(int i=0;i<colnum_objects;i+=2)
		// 	{
			
		// 		int te=*((int*)(sorted_cols+i+1));
		// 		cout<<pid<<": "<<cr_data[i]<<"  "<<cr_data[i+1]<<"->  "<<sorted_cols[i]<<"  "<<te<< endl;
		// 	}	
		// }


		if(pid==0)
		{
			for(int i=0;i<2*num_entries;i+=2)
			{
				all_nodes[i/2].valfloat=sorted_cols[i];
				int te= *((int*)(sorted_cols+i+1));
				all_nodes[i/2].valint= te;
			}
			sort(all_nodes.begin(),all_nodes.end(), structrowcomparator);	
		}


		free(sorted_cols);
		colelem_count.clear();
		free(colelem_count_prefix);
		free(sorted_rows);
		free(colrecv_data);
		free(colall_data);
		free(cr_data);
		free(recv_cr);
		// free(all_data);
		// free(rc_data);
		free(recv_data);
		free(recv_rc);
		free(elem_count_prefix);
		elem_count.clear();

		// cout<<"a looper is finished"<<endl;

		// exit(1);

	}

	// final writing
	if(pid==0)
	{

		FILE *file = fopen(argv[2],"wb") ;
		for(int i=0 ; i<num_entries ; i++)
		{
			
			fwrite(&all_nodes[i].r,4,1,file) ;
			fwrite(&all_nodes[i].c,4,1,file) ;
			fwrite(&all_nodes[i].valint,4,1,file) ;
			fwrite(&all_nodes[i].valfloat,4,1,file) ;
			
		}
		fclose(file);
	}

	// cout<<"file writing done"<<endl;
	


MPI_Finalize();
	return 0;
	}

	
