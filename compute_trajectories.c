#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

//allocate memory for new matrix
float **new_matrix(int N, int M, float length, int seed);
//free memory of matrix L
void free_matrix(float **matrix, int N);
//print Matrix L using + for +1 and - for -1
void print_matrix(float **matrix, int N, int M);
//perform sim_length integration step with VV algorithm
void integrate_eom(float **phase_point, int sim_length, float length, int N, float dt, float m, float G, float radius, float G_cutoff);
//check if two balls overlap and update their velocities accordingly
void collision_check(float **phase_point, int N, float radius);


int main(){
    
    int sim_length = 3000;
    int N = 500;
    float dt = .01;
    float G = 0.05;
    float m = 1;
    float length = 10;
    float radius = 0.1;
    float r_cutoff = 10;
    float G_cutoff = radius;
    
    int seed = 1;


    float **phase_point = new_matrix(N,4,length, seed);


    printf("Computing... \n");
    integrate_eom(phase_point,sim_length,length,N,dt,m,G,radius,G_cutoff);
    printf("Done\n");
    free_matrix(phase_point, N);

    return 0;
}



void integrate_eom(float **phase_point, int sim_length, float length, int N, float dt, float m, float G, float radius, float G_cutoff){
    
    int i,j,t=0;
    float F_x_temp[N];
    float F_y_temp[N];

    //create an array of pointers to files; one for each particle
    FILE** output_files = malloc(sizeof(FILE*) * N);

    
    //create filenames that include the particle index i and open them
    for (i=0;i<N;i++){
        char* data = "data\\data";
        char* extension = ".txt";

        char filename[strlen(data)+strlen(extension)+1];


        snprintf( filename, sizeof( filename ) + sizeof( "_" ) + sizeof(i), "%s_%d%s", data,i, extension );

        output_files[i] = fopen(filename,"w");
    }
    
    int progress;
    int n_print = (int) sim_length/100;
    while (t < sim_length){

        if (t % n_print == 0){
            progress =(int) (double)t/(double)sim_length * 100;
            
            printf("%d%%\r",progress);
        }

        collision_check(phase_point, N, radius);

        //First, update ALL positions
        for (i=0;i<N;i++){

            float x = phase_point[i][0];
            float y = phase_point[i][1];
            float v_x = phase_point[i][2];
            float v_y = phase_point[i][3];

            float F_x = 0;
            float F_y = 0;

            //sum up all the forces from particles j on particle i
            for (j=0; j<N;j++){
                //no self interaction
                if (i != j){
                    float d_x = x - phase_point[j][0];
                    float d_y = y - phase_point[j][1];

                    //extend force accors pbc
                    if (d_x > length){
                        d_x -= 2*length;
                    }
                    else if (d_x < -length){
                        d_x += 2*length;
                    }
                    if (d_y > length){
                        d_y -= 2*length;
                    }
                    else if (d_y < -length){
                        d_y += 2*length;
                    }

                    float r = sqrt(d_x*d_x + d_y*d_y);
                    if (r > G_cutoff){
                        F_x += -G*m*m/(r*r*r)*d_x;
                        F_y += -G*m*m/(r*r*r)*d_y;
                    }
                }
            }

            //save current forces on particle i for calculating velocities later
            F_x_temp[i]=F_x;
            F_y_temp[i]=F_y;

            float x_new = x + v_x*dt + F_x/(2*m)*dt*dt;
            float y_new = y + v_y*dt + F_y/(2*m)*dt*dt;
            
            phase_point[i][0] = x_new;
            phase_point[i][1] = y_new;
        }


        //Then, update the velocities
        for (i=0;i<N;i++){
        
            float x = phase_point[i][0];
            float y = phase_point[i][1];
            float v_x = phase_point[i][2];
            float v_y = phase_point[i][3];

            float F_x = 0;
            float F_y = 0;

            float F_x_new = F_x_temp[i];
            float F_y_new = F_y_temp[i];

            //sum up all the forces from particles j on particle i
            for (j=0; j<N;j++){
                //no self interaction
                if (i != j){
                    float d_x = x - phase_point[j][0];
                    float d_y = y - phase_point[j][1];

                    //extend force accors pbc
                    if (d_x > length){
                        d_x -= 2*length;
                    }
                    else if (d_x < -length){
                        d_x += 2*length;
                    }
                    if (d_y > length){
                        d_y -= 2*length;
                    }
                    else if (d_y < -length){
                        d_y += 2*length;
                    }

                    float r = sqrt(d_x*d_x + d_y*d_y);

                    if (r > G_cutoff){
                        F_x += -G*m*m/(r*r*r)*d_x;
                        F_y += -G*m*m/(r*r*r)*d_y;
                    }
                }
            }



            float v_x_new = v_x + (F_x+F_x_new)/(2*m)*dt; 
            float v_y_new = v_y + (F_y+F_y_new)/(2*m)*dt;

            phase_point[i][2] = v_x_new;
            phase_point[i][3] = v_y_new;
            

            //periodic boundary conditions
            if (phase_point[i][0] > length){
                phase_point[i][0] -= 2*length;
            }
            else if (phase_point[i][0] < (-length)){
                phase_point[i][0] += 2*length;
            }

            if (phase_point[i][1] > length){
                phase_point[i][1] -= 2*length;
            }
            else if (phase_point[i][1] < (-length)){
                phase_point[i][1] += 2*length;
            }

        }
        //write the current configurations to the text files
        for (i=0;i<N;i++){
            fprintf(output_files[i],"%f %f\n",phase_point[i][0],phase_point[i][1]);
        }
        
        t++;
    }

    //close and free the files
    for (i=0; i<N;i++){
        fclose(output_files[i]);
    }

    free(output_files);
}





void collision_check(float **phase_point, int N, float radius){
    int i,j;

    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            float d_x = phase_point[i][0]-phase_point[j][0];
            float d_y = phase_point[i][1]-phase_point[j][1];

            float d_abs_squared = d_x*d_x+d_y*d_y;

            //check whether ball i and ball j are overlapping
            if (d_abs_squared < (4*radius*radius)){
                float delta_v_x = phase_point[j][2]-phase_point[i][2];
                float delta_v_y = phase_point[j][3]-phase_point[i][3];

                float delta_v_dot = d_x*delta_v_x + d_y*delta_v_y;

                //check whether ball i and ball j are moving towards each other (otherwise they get stuck)
                if (delta_v_dot > 0){
                    //distance between ball i and j
                    float d_abs = sqrt(d_abs_squared);
                    d_x = d_x/d_abs;
                    d_y = d_y/d_abs;
                    float a = d_x*(phase_point[i][2]-phase_point[j][2]) + d_y*(phase_point[i][3]-phase_point[j][3]);
                
                    phase_point[i][2] -= a*d_x;
                    phase_point[i][3] -= a*d_y;

                    phase_point[j][2] += a*d_x;
                    phase_point[j][3] += a*d_y;

                    phase_point[j][0] = phase_point[i][0] - 2*radius*d_x;
                    phase_point[j][1] = phase_point[i][1] - 2*radius*d_y;
                }
            }
        }
    }
}





//allocate memory for new matrix; each row represents a particle with columns [x,y,vx,vy]
float **new_matrix(int N, int M, float length, int seed){

	int i;
	float r;
	float **matrix = (float **)malloc(N * sizeof(float *));
	
	srand(seed);

	for (i=0; i<N; i++){
		matrix[i] = (float *)malloc(M * sizeof(float *));
	}


	//set up random initial configuration with zero velocities
	for (i=0; i<N; i++){

        r = (float) rand() / (double) RAND_MAX;
		matrix[i][0] = length*(2*r-1);

        r = (float) rand() / (double) RAND_MAX;
        matrix[i][1] = length*(2*r-1);

        matrix[i][2] = 0;
        matrix[i][3] = 0;

	}
	return matrix;
}


//free memory that was allocated for matrix L
void free_matrix(float **matrix, int N){

	int i;

	for (i=0; i<N; i++){
		free(matrix[i]);
	}
}


//print matrix
void print_matrix(float **matrix, int N, int M){

	int i,j;

	for (i=0; i<N; i++){
		for (j=0; j<M; j++){

			printf("%f ",matrix[i][j]);
		}
		printf("\n");
	}
}
