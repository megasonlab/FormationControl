//Parameter.h for 2D Vertex Model
#ifndef _PARAMETER_H
#define _PARAMETER_H

#define Pi (3.14159265358979323846)

//GLSC3D Setting
#define WINDOW_SIZE_X (2000)
#define WINDOW_SIZE_Y (1200)

#define usecapture (0)
#define GLSC_INTV (1000)

#define Max_VertexNumber (100000)
#define Max_CellNumber (1000)

#define Max_Time_Step (100000)

//Time discritization width
double dt = 1.0e-3;
//Time
double t = 0.0;

//Drawing width
double L = 250;


//When the variable "Cell_Adhesion" is defined,
//posterior FP and HC cells are connected with posterior notochord cells
#define Cell_Adhesion (0)
int Divide_Step = 800;
//int Divide_Step = 1200;


int Init_Cell_Number = 32;

int FP_Number = Init_Cell_Number;
int HC_Number = Init_Cell_Number;
int NOTO_Number = Init_Cell_Number*2;

int FP_Vertex_Number = 2 * (Init_Cell_Number + 1);
int HC_Vertex_Number = 2 * (Init_Cell_Number + 1);
int NOTO_Vertex_Number = 2 * (Init_Cell_Number*2 + 1);

double Cell_Width = 2.0;
double NOTO_Cell_Width = 1.0;

double NOTO_Height = 2.0;
double FP_Height = 1.0;
double HC_Height = 1.0;


//Volume Conservation
double a = 10.0;
double Ideal_NOTO_Area = 2.0;
double Ideal_FP_Area = 2.0;
double Ideal_HC_Area = 2.0;

//Edge Length Conservation
double b = 0.5;
double Init_Length = 2.0;
double Ideal_Length = 0.0;

//Adhesion to ECM: Parameter c in Table 1
double Basal_Adhesion_Parameter = 100.0;

//Repulsion for Apical surface: Parameter d in Table 1
double Apical_Repulsion_Parameter = 100.0;

//Parameter sigma in Table 1
double sigma = 0.2;

//Potential energy threshold: Parameter Uth in Table 1
double Threshold = 2.5;

//Migration Force Parameter
double Velocity_Index = 1.0; //Parameter alpha
double Migration_Parameter = 2.0; //Parameter v_max in Table 1

//Cell Division Timing: Parameter S_i in Table 1
double Divide_Timing = 4.0;

//Time constant for a differential equation of Migration Force
double eta = 1.0;

//When the following time is elapsed, division cells are colored
double Coloring_Time = Max_Time_Step * dt * 0.5;

//Parameter for controlling x-direction movement for the left side cells
double Control_Parameter = 1.0e-15;

#endif
