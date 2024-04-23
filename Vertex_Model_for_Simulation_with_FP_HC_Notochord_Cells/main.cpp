#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

//GLSC3D for Drawing Simulation Result
//https://github.com/GLSC3DProject/GLSC3D
#include<glsc3d_3.h>

//Mersenne Twister
//http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html
extern "C" {
    #include <mt19937ar.h>
};

//Header files for Simulation of Vertex model
#include"./Header/Vector.h"
#include"./Header/Parameter.h"
#include"./Header/Function.h"

//Main Function
int main(int argc, char *argv[]) {
    
    //Command Line Parameters
    Velocity_Index = atof(argv[1]);
    Migration_Parameter = atof(argv[2]);
    int test_index = atoi(argv[3]);
    //Example: We set Velocity_Index = 1, Migration_Parameter = 2, and test_index = 1.
    //We execute "makefile", then an executable file "run" is generated.
    //Next we excute the executable File "run" as "./run 1 2 1".
    
    //Seed of Random Number
//    unsigned int seed = (unsigned int)time(NULL) + test_index;
    unsigned int seed = 1;
    init_genrand(seed);
    
    FILE *fp;
    char filename[256];

    //Array of Vertex
    Vector2D *FP_Vertex = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    Vector2D *HC_Vertex = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    Vector2D *NOTO_Vertex = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    
    Vector2D *New_FP_Vertex = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    Vector2D *New_HC_Vertex = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    Vector2D *New_NOTO_Vertex = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    
    Vector2D *FP_Force = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    Vector2D *HC_Force = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    Vector2D *NOTO_Force = (Vector2D *)malloc(sizeof(Vector2D) * Max_VertexNumber);
    
    //Initialize of Vertices
    for(int i=0; i<=Init_Cell_Number; i++) {
        FP_Vertex[i+0*(Init_Cell_Number+1)]   = {i * Cell_Width, NOTO_Height * 0.5};
        FP_Vertex[i+1*(Init_Cell_Number+1)]   = {i * Cell_Width, NOTO_Height * 0.5 + FP_Height};
        HC_Vertex[i+0*(Init_Cell_Number+1)]   = {i * Cell_Width, -NOTO_Height * 0.5};
        HC_Vertex[i+1*(Init_Cell_Number+1)]   = {i * Cell_Width, -NOTO_Height * 0.5 - HC_Height};
    }
    for(int i=0; i<=Init_Cell_Number * 2; i++) {
        NOTO_Vertex[i+0*(Init_Cell_Number*2+1)] = {i * NOTO_Cell_Width, -NOTO_Height * 0.5};
        NOTO_Vertex[i+1*(Init_Cell_Number*2+1)] = {i * NOTO_Cell_Width, NOTO_Height * 0.5};
    }

    //Indices of Fixed Vertices
    int FixedVertexBottom = 0;
    int FixedVertexTop = 0+1*(Init_Cell_Number+1);
    int FixedVertexBottom_NOTO = 0;
    int FixedVertexTop_NOTO = 0+1*(Init_Cell_Number*2+1);
    
    //FP Cell Informations
    Cell *FP_Cell = (Cell *)malloc(sizeof(Cell) * Max_VertexNumber);
    for(int i=0; i<Init_Cell_Number; i++) {
        FP_Cell[i].number = i;
        FP_Cell[i].index[0] = i;
        FP_Cell[i].index[1] = i+1;
        FP_Cell[i].index[2] = i+1+1*(Init_Cell_Number+1);
        FP_Cell[i].index[3] = i+1*(Init_Cell_Number+1);
        FP_Cell[i].flag = 0;
        FP_Cell[i].divide_flag = 0;
        FP_Cell[i].time = 0;
    }
    //HC Cell Informations
    Cell *HC_Cell = (Cell *)malloc(sizeof(Cell) * Max_VertexNumber);
    for(int i=0; i<Init_Cell_Number; i++) {
        HC_Cell[i].number = i;
        HC_Cell[i].index[3] = i;
        HC_Cell[i].index[2] = i+1;
        HC_Cell[i].index[1] = i+1+1*(Init_Cell_Number+1);
        HC_Cell[i].index[0] = i+1*(Init_Cell_Number+1);
        HC_Cell[i].flag = 0;
        HC_Cell[i].divide_flag = 0;
        HC_Cell[i].time = 0;
    }
    //NOTO Cell Informations
    Cell *NOTO_Cell = (Cell *)malloc(sizeof(Cell) * Max_VertexNumber);
    for(int i=0; i<Init_Cell_Number*2; i++) {
        NOTO_Cell[i].number = i;
        NOTO_Cell[i].index[0] = i;
        NOTO_Cell[i].index[1] = i+1;
        NOTO_Cell[i].index[2] = i+1+1*(Init_Cell_Number*2+1);
        NOTO_Cell[i].index[3] = i+1*(Init_Cell_Number*2+1);
        NOTO_Cell[i].flag = 0;
        NOTO_Cell[i].divide_flag = 0;
        NOTO_Cell[i].time = 0;
    }

    //Potential Energy of Apical and Basal Edges
    double *FP_Cell_Energy_Data = (double *)malloc(sizeof(double) * Max_CellNumber);
    double *HC_Cell_Energy_Data = (double *)malloc(sizeof(double) * Max_CellNumber);
    for(int i=0; i<Max_CellNumber; i++) { FP_Cell_Energy_Data[i] = 0.0; }
    for(int i=0; i<Max_CellNumber; i++) { HC_Cell_Energy_Data[i] = 0.0; }
    
    //Array for Cell Label
    int *FP_Cell_Order = (int *)malloc(sizeof(int) * Max_CellNumber);
    int *tmp_FP_Cell_Order = (int *)malloc(sizeof(int) * Max_CellNumber);
    for(int i=0; i<FP_Number; i++) { FP_Cell_Order[i] = i; }

    int *HC_Cell_Order = (int *)malloc(sizeof(int) * Max_CellNumber);
    int *tmp_HC_Cell_Order = (int *)malloc(sizeof(int) * Max_CellNumber);
    for(int i=0; i<HC_Number; i++) { HC_Cell_Order[i] = i; }
    
    //Cell Centroid
    Vector2D *FP_CellCentroid = (Vector2D *)malloc(sizeof(Vector2D) * Max_CellNumber);
    Vector2D *HC_CellCentroid = (Vector2D *)malloc(sizeof(Vector2D) * Max_CellNumber);
    
    //Array for Migration Force
    double *FP_Cell_Mig_Force = (double *)malloc(sizeof(double) * Max_CellNumber);
    double *HC_Cell_Mig_Force = (double *)malloc(sizeof(double) * Max_CellNumber);
    for(int i=0; i<FP_Number; i++) { FP_Cell_Mig_Force[i] = 0.0; }
    for(int i=0; i<HC_Number; i++) { HC_Cell_Mig_Force[i] = 0.0; }

    //GLSC Setting
    g_enable_highdpi();
    g_set_antialiasing(4);
    
    //Window Setting
//    g_init("Simulation of 2D Vertex model", WINDOW_SIZE_X, WINDOW_SIZE_Y);
    g_init(G_OFF_SCREEN, WINDOW_SIZE_X, WINDOW_SIZE_Y);
    
    //Background Color
    g_scr_color(1,1,1);
    
    //GLSC Screen Setting
    g_def_scale_2D(0, -1.0, L, -10, 10, 0, 0, WINDOW_SIZE_X, WINDOW_SIZE_Y/3);
    g_def_scale_2D(1, -1.0, L, -10, 10, 0, WINDOW_SIZE_Y/3, WINDOW_SIZE_X, WINDOW_SIZE_Y/3);
    g_def_scale_2D(2, -1.0, L, -1.0, 2.0 * Threshold, 0, WINDOW_SIZE_Y/3*2, WINDOW_SIZE_X, WINDOW_SIZE_Y/3);

#ifdef usecapture
    sprintf(filename, "mkdir Image_%d", test_index);
    system(filename);
    sprintf(filename, "./Image_%d/Capture_%4.4f_%4.4f", test_index, Velocity_Index, Migration_Parameter);
    g_capture_set(filename);
#endif
    
    sprintf(filename, "mkdir Data_%d", test_index);
    system(filename);
    //FP Cell Division Data
    sprintf(filename, "./Data_%d/FP_Division_Data_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
    fp = fopen(filename, "a");
    fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
    fclose(fp);

    //HC Cell Division Data
//    sprintf(filename, "./Data_%d/HC_Division_Data_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
//    fp = fopen(filename, "a");
//    fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
//    fclose(fp);

    sprintf(filename, "mkdir Angle_Data_%d", test_index);
    system(filename);
    
    sprintf(filename, "mkdir Migration_Activity_Data_%d", test_index);
    system(filename);
    
    sprintf(filename, "./Data_%d/Seed_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
    fp = fopen(filename, "a");
    fprintf(fp, "%u\n", seed);
    fclose(fp);
    
    //Main Time Loop of Numerical Simulation
    for (int i_time = 0; i_time <= Max_Time_Step; i_time++) {
        
        //Current Time
        t = i_time * dt;
        
        //Drawing by GLSC3D
        if (i_time % GLSC_INTV == 0) {
            g_cls();

            //Result of FP Cells
            g_sel_scale(0);
            g_line_width(0.5);
            g_line_color(0,0,0,1);
            g_boundary();

            //Division FP Cells are drawn by Green Color
            for (int i=0; i<FP_Number; i++) {
                g_area_color(0.8, 0.8, 0.8, 1);
                if (FP_Cell[i].flag == 1) g_area_color(0.0, 0.8, 0.0, 0.5);
                
                for (int j=0; j<4; j++) {
                    Vector2D tmp_Centroid = Zero2D;
                    for (int j=0; j<4; j++) { tmp_Centroid += FP_Vertex[ FP_Cell[i].index[j] ] * 0.25; }
                    g_triangle_2D(FP_Vertex[FP_Cell[i].index[j]].x, FP_Vertex[FP_Cell[i].index[j]].y,
                                  FP_Vertex[FP_Cell[i].index[(j+1)%4]].x, FP_Vertex[FP_Cell[i].index[(j+1)%4]].y,
                                  tmp_Centroid.x, tmp_Centroid.y, G_NO, G_YES);
                }
            }
            
            for (int i=0; i<FP_Number; i++) {
                g_line_color(0,0,0,1);
                g_line_width(1.0);
                g_move_2D(FP_Vertex[FP_Cell[i].index[0]].x, FP_Vertex[FP_Cell[i].index[0]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[3]].x, FP_Vertex[FP_Cell[i].index[3]].y);
                g_move_2D(FP_Vertex[FP_Cell[i].index[2]].x, FP_Vertex[FP_Cell[i].index[2]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[1]].x, FP_Vertex[FP_Cell[i].index[1]].y);
                
                double tmp_value = (FP_Cell_Energy_Data[i]) / (2.0 * Threshold);
                g_line_color(tmp_value, 0.0, 1.0 - tmp_value, 1);
                g_line_width(2.0);
                g_move_2D(FP_Vertex[FP_Cell[i].index[0]].x, FP_Vertex[FP_Cell[i].index[0]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[1]].x, FP_Vertex[FP_Cell[i].index[1]].y);
                g_move_2D(FP_Vertex[FP_Cell[i].index[2]].x, FP_Vertex[FP_Cell[i].index[2]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[3]].x, FP_Vertex[FP_Cell[i].index[3]].y);
            }
            
            //Result with NOTO and HC Cells
            g_sel_scale(1);
            g_line_width(0.5);
            g_line_color(0,0,0,1);
            g_boundary();

            //Division FP Cells are drawn by Green Color
            for (int i=0; i<FP_Number; i++) {
                g_area_color(0.8, 0.8, 0.8, 1);
                if (FP_Cell[i].flag == 1) g_area_color(0.0, 0.8, 0.0, 0.5);
                Vector2D tmp_Centroid = Zero2D;
                for (int j=0; j<4; j++) { tmp_Centroid += FP_Vertex[ FP_Cell[i].index[j] ] * 0.25; }
                for (int j=0; j<4; j++) {
                    g_triangle_2D(FP_Vertex[FP_Cell[i].index[j]].x, FP_Vertex[FP_Cell[i].index[j]].y,
                                  FP_Vertex[FP_Cell[i].index[(j+1)%4]].x, FP_Vertex[FP_Cell[i].index[(j+1)%4]].y,
                                  tmp_Centroid.x, tmp_Centroid.y, G_NO, G_YES);
                }
            }
            for (int i=0; i<FP_Number; i++) {
                g_line_color(0,0,0,1);
                g_line_width(1.0);
                g_move_2D(FP_Vertex[FP_Cell[i].index[0]].x, FP_Vertex[FP_Cell[i].index[0]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[3]].x, FP_Vertex[FP_Cell[i].index[3]].y);
                g_move_2D(FP_Vertex[FP_Cell[i].index[2]].x, FP_Vertex[FP_Cell[i].index[2]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[1]].x, FP_Vertex[FP_Cell[i].index[1]].y);
                
                double tmp_value = (FP_Cell_Energy_Data[i]) / (2.0 * Threshold);
                g_line_color(tmp_value, 0.0, 1.0 - tmp_value, 1);
                g_line_width(2.0);
                g_move_2D(FP_Vertex[FP_Cell[i].index[0]].x, FP_Vertex[FP_Cell[i].index[0]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[1]].x, FP_Vertex[FP_Cell[i].index[1]].y);
                g_move_2D(FP_Vertex[FP_Cell[i].index[2]].x, FP_Vertex[FP_Cell[i].index[2]].y);
                g_plot_2D(FP_Vertex[FP_Cell[i].index[3]].x, FP_Vertex[FP_Cell[i].index[3]].y);
            }
            
            //Division HC Cells are drawn by Green Color
            for (int i=0; i<HC_Number; i++) {
                g_area_color(0.8, 0.8, 0.8, 1);
                if (HC_Cell[i].flag == 1) g_area_color(0.0, 0.8, 0.0, 0.5);
                Vector2D tmp_Centroid = Zero2D;
                for (int j=0; j<4; j++) { tmp_Centroid += HC_Vertex[ HC_Cell[i].index[j] ] * 0.25; }
                for (int j=0; j<4; j++) {
                    g_triangle_2D(HC_Vertex[HC_Cell[i].index[j]].x, HC_Vertex[HC_Cell[i].index[j]].y,
                                  HC_Vertex[HC_Cell[i].index[(j+1)%4]].x, HC_Vertex[HC_Cell[i].index[(j+1)%4]].y,
                                  tmp_Centroid.x, tmp_Centroid.y, G_NO, G_YES);
                }
            }
            for (int i=0; i<HC_Number; i++) {
                g_line_color(0,0,0,1);
                g_line_width(1.0);
                g_move_2D(HC_Vertex[HC_Cell[i].index[0]].x, HC_Vertex[HC_Cell[i].index[0]].y);
                g_plot_2D(HC_Vertex[HC_Cell[i].index[3]].x, HC_Vertex[HC_Cell[i].index[3]].y);
                g_move_2D(HC_Vertex[HC_Cell[i].index[2]].x, HC_Vertex[HC_Cell[i].index[2]].y);
                g_plot_2D(HC_Vertex[HC_Cell[i].index[1]].x, HC_Vertex[HC_Cell[i].index[1]].y);
                
                double tmp_value = (FP_Cell_Energy_Data[i]) / (2.0 * Threshold);
                g_line_color(tmp_value, 0.0, 1.0 - tmp_value, 1);
                g_line_width(2.0);
                g_move_2D(HC_Vertex[HC_Cell[i].index[0]].x, HC_Vertex[HC_Cell[i].index[0]].y);
                g_plot_2D(HC_Vertex[HC_Cell[i].index[1]].x, HC_Vertex[HC_Cell[i].index[1]].y);
                g_move_2D(HC_Vertex[HC_Cell[i].index[2]].x, HC_Vertex[HC_Cell[i].index[2]].y);
                g_plot_2D(HC_Vertex[HC_Cell[i].index[3]].x, HC_Vertex[HC_Cell[i].index[3]].y);
            }
            
            //NOTO Cells are drawn by Gray Color
            for (int i=0; i<NOTO_Number; i++) {
                g_area_color(0.8, 0.8, 0.8, 1);
                Vector2D tmp_Centroid = Zero2D;
                for (int j=0; j<4; j++) { tmp_Centroid += NOTO_Vertex[ NOTO_Cell[i].index[j] ] * 0.25; }
                for (int j=0; j<4; j++) {
                    g_triangle_2D(NOTO_Vertex[NOTO_Cell[i].index[j]].x, NOTO_Vertex[NOTO_Cell[i].index[j]].y,
                                  NOTO_Vertex[NOTO_Cell[i].index[(j+1)%4]].x, NOTO_Vertex[NOTO_Cell[i].index[(j+1)%4]].y,
                                  tmp_Centroid.x, tmp_Centroid.y, G_NO, G_YES);
                }
            }
            for (int i=0; i<NOTO_Number; i++) {
                g_line_color(0,0,0,1);
                g_line_width(1.0);
                g_move_2D(NOTO_Vertex[NOTO_Cell[i].index[0]].x, NOTO_Vertex[NOTO_Cell[i].index[0]].y);
                g_plot_2D(NOTO_Vertex[NOTO_Cell[i].index[3]].x, NOTO_Vertex[NOTO_Cell[i].index[3]].y);
                g_move_2D(NOTO_Vertex[NOTO_Cell[i].index[2]].x, NOTO_Vertex[NOTO_Cell[i].index[2]].y);
                g_plot_2D(NOTO_Vertex[NOTO_Cell[i].index[1]].x, NOTO_Vertex[NOTO_Cell[i].index[1]].y);
                g_move_2D(NOTO_Vertex[NOTO_Cell[i].index[0]].x, NOTO_Vertex[NOTO_Cell[i].index[0]].y);
                g_plot_2D(NOTO_Vertex[NOTO_Cell[i].index[1]].x, NOTO_Vertex[NOTO_Cell[i].index[1]].y);
                g_move_2D(NOTO_Vertex[NOTO_Cell[i].index[2]].x, NOTO_Vertex[NOTO_Cell[i].index[2]].y);
                g_plot_2D(NOTO_Vertex[NOTO_Cell[i].index[3]].x, NOTO_Vertex[NOTO_Cell[i].index[3]].y);
            }
            
            //Angle Data at ECM Side Posterior Vertex of a Sinle FP Cell
            sprintf(filename, "./Angle_Data_%d/FP_Angle_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
            fp = fopen(filename, "w");
            fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
            for (int i=0; i<FP_Number; i++) {
                Vector2D tmp_next = FP_Vertex[FP_Cell[i].index[2]] - FP_Vertex[FP_Cell[i].index[1]];
                tmp_next = tmp_next / ~tmp_next;
                Vector2D tmp_prev = FP_Vertex[FP_Cell[i].index[0]] - FP_Vertex[FP_Cell[i].index[1]];
                tmp_prev = tmp_prev / ~tmp_prev;
                double tmp_angle = tmp_prev*tmp_next;
                fprintf(fp, "%f %f\n", FP_Vertex[FP_Cell[i].index[1]].x, tmp_angle);
            }
           fclose(fp);

            //Angle Data at ECM Side Posterior Vertex of a Sinle HC Cell
//           sprintf(filename, "./Angle_Data_%d/HC_Angle_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
//            fp = fopen(filename, "w");
//            fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
//            for (int i=0; i<HC_Number; i++) {
//                Vector2D tmp_next = HC_Vertex[HC_Cell[i].index[3]] - HC_Vertex[HC_Cell[i].index[2]];
//                tmp_next = tmp_next / ~tmp_next;
//                Vector2D tmp_prev = HC_Vertex[HC_Cell[i].index[1]] - HC_Vertex[HC_Cell[i].index[2]];
//                tmp_prev = tmp_prev / ~tmp_prev;
//                double tmp_angle = tmp_prev*tmp_next;
//                fprintf(fp, "%f %f\n", HC_Vertex[HC_Cell[i].index[1]].x, tmp_angle);
//            }
//           fclose(fp);
            
            //Migration Activity Data of a Sinle FP Cell
//            sprintf(filename, "./Migration_Activity_Data_%d/Migration_Activity_%4.4f_%4.4f_%4.4f.dat", t, test_index, Velocity_Index, Migration_Parameter);
//            fp = fopen(filename, "w");
//            fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
//            for (int i=0; i<FP_Number; i++) {
//               fprintf(fp, "%f %f\n", FP_Vertex[FP_Cell[i].index[1]].x,
//                       FP_Cell_Mig_Force_w[i]);
//            }
//           fclose(fp);
            
            //Result of Graphs
            g_sel_scale(2);
            g_line_width(1.0);
            g_line_color(0,0,0,1);
            g_boundary();
            
            g_text_size(10);
            g_text_2D_virtual(15.0, 2.0 * Threshold * 0.8, "t = %2.2f", t);
            g_text_2D_virtual(25.0, 2.0 * Threshold * 0.8, "FP_Cell_Number = %d", FP_Number);
            g_text_2D_virtual(45.0, 2.0 * Threshold * 0.8, "HC_Cell_Number = %d", HC_Number);
            g_text_2D_virtual(65.0, 2.0 * Threshold * 0.8, "NOTO_Cell_Number = %d", NOTO_Number);
            
            g_line_color(0.5, 0.5, 0.5, 1.0);
            g_move_2D(-10.0, 0.0); g_plot_2D(L, 0.0);
            g_move_2D(0.0, -10.0); g_plot_2D(0.0, 2.0*Threshold);
            g_line_type(1);
            g_line_color(0.0, 0.0, 0.5, 1.0);
            g_move_2D(-10.0, Threshold); g_plot_2D(L, Threshold);
            g_line_type(0);
            
            g_marker_radius(0.05);
            g_marker_type(1);
            g_marker_color(0.5, 0.5, 0.5, 1.0);
            for (int i=0; i<FP_Number; i++) {
                g_marker_2D(FP_CellCentroid[FP_Cell_Order[i]].x, FP_Cell_Energy_Data[FP_Cell_Order[i]]);
            }

            g_line_color(0.0, 0.0, 1.0, 1.0);
            g_line_width(2.0);
            for (int i=0; i<FP_Number-1; i++) {
                g_move_2D(FP_CellCentroid[FP_Cell_Order[i]].x, FP_Cell_Energy_Data[FP_Cell_Order[i]]);
                g_plot_2D(FP_CellCentroid[FP_Cell_Order[i+1]].x, FP_Cell_Energy_Data[FP_Cell_Order[i+1]]);
            }
            
            g_marker_radius(0.05);
            g_marker_type(1);
            g_marker_color(0.5, 0.5, 0.5, 1.0);
            for (int i=0; i<FP_Number; i++) {
                g_marker_2D(FP_CellCentroid[FP_Cell_Order[i]].x, FP_Cell_Mig_Force[FP_Cell_Order[i]]);
            }
            
            g_line_color(1.0, 0.0, 1.0, 1.0);
            g_line_width(2.0);
            for (int i=0; i<FP_Number-1; i++) {
                g_move_2D(FP_CellCentroid[FP_Cell_Order[i]].x, FP_Cell_Mig_Force[FP_Cell_Order[i]]);
                g_plot_2D(FP_CellCentroid[FP_Cell_Order[i+1]].x, FP_Cell_Mig_Force[FP_Cell_Order[i+1]]);
            }
        
            g_finish();
            
#ifdef usecapture
            g_capture();
#endif
        }
        
        //If the potential energy for edge length (Apical and Basal) is greater than a specified threshold,
        //a time to cell division is calculated.
        //FP Cell Calculation
        for (int i=0; i<FP_Number; i++)
            if (FP_Cell_Energy_Data[FP_Cell_Order[i]] > Threshold)
                FP_Cell[FP_Cell_Order[i]].divide_flag = 1;
        for (int i=0; i<FP_Number; i++)
            if (FP_Cell[FP_Cell_Order[i]].divide_flag == 1)
                FP_Cell[FP_Cell_Order[i]].time += dt;
        //HC Cell Calculation
        for (int i=0; i<HC_Number; i++)
            if (HC_Cell_Energy_Data[HC_Cell_Order[i]] > Threshold)
                HC_Cell[HC_Cell_Order[i]].time += dt;
        for (int i=0; i<HC_Number; i++)
            if (HC_Cell[HC_Cell_Order[i]].divide_flag == 1)
                HC_Cell[HC_Cell_Order[i]].time += dt;
        
        double Max_FP_Centroid_x = FP_CellCentroid[0].x;
        double Min_FP_Centroid_x = FP_CellCentroid[0].x;
        for (int i=1; i<FP_Number; i++) {
            if (Max_FP_Centroid_x < FP_CellCentroid[i].x) Max_FP_Centroid_x = FP_CellCentroid[i].x;
            if (Min_FP_Centroid_x > FP_CellCentroid[i].x) Min_FP_Centroid_x = FP_CellCentroid[i].x;
        }
        double Max_HC_Centroid_x = HC_CellCentroid[0].x;
        double Min_HC_Centroid_x = HC_CellCentroid[0].x;
        for (int i=1; i<HC_Number; i++) {
            if (Max_HC_Centroid_x < HC_CellCentroid[i].x) Max_HC_Centroid_x = HC_CellCentroid[i].x;
            if (Min_HC_Centroid_x > HC_CellCentroid[i].x) Min_HC_Centroid_x = HC_CellCentroid[i].x;
        }
        
        //FP Cell Division
        for(int i=0; i<FP_Number; i++) {
            tmp_FP_Cell_Order[i] = FP_Cell_Order[i];
        }
        int tmp_CellNumber = FP_Number;
        int CellNumber_Counter = 0;
        for (int i=0; i<FP_Number; i++) {
            if ((FP_Cell_Energy_Data[tmp_FP_Cell_Order[i]] > Threshold &&
                 FP_Cell[tmp_FP_Cell_Order[i]].time > Divide_Timing)
                && FP_Number < Max_CellNumber) {
                
                double tmp_FP_Centroid_x = FP_CellCentroid[tmp_FP_Cell_Order[i]].x;
                double tmp_length = (tmp_FP_Centroid_x - Min_FP_Centroid_x) / (Max_FP_Centroid_x - Min_FP_Centroid_x);

                Vector2D Bottom_Vertex = (FP_Vertex[ FP_Cell[tmp_FP_Cell_Order[i]].index[0] ]
                                          + FP_Vertex[ FP_Cell[tmp_FP_Cell_Order[i]].index[1] ]) * 0.5;
                Vector2D Top_Vertex = (FP_Vertex[ FP_Cell[tmp_FP_Cell_Order[i]].index[2] ]
                                       + FP_Vertex[ FP_Cell[tmp_FP_Cell_Order[i]].index[3] ]) * 0.5;
                FP_Vertex[FP_Vertex_Number++] = Bottom_Vertex;
                FP_Vertex[FP_Vertex_Number++] = Top_Vertex;

                int index[4];
                for (int j=0; j<4; j++) index[j] = FP_Cell[tmp_FP_Cell_Order[i]].index[j];
                
                //Left Cell
                FP_Cell[FP_Number].number = 4;
                FP_Cell[FP_Number].index[0] = index[0];
                FP_Cell[FP_Number].index[1] = FP_Vertex_Number-2;
                FP_Cell[FP_Number].index[2] = FP_Vertex_Number-1;
                FP_Cell[FP_Number].index[3] = index[3];
                FP_Cell_Mig_Force[FP_Number] = 0.0;

                if (t > Coloring_Time) FP_Cell[tmp_FP_Cell_Order[i]].flag = 1;
                else FP_Cell[tmp_FP_Cell_Order[i]].flag = 0;
                FP_Cell[tmp_FP_Cell_Order[i]].divide_flag = 0;
                FP_Cell[tmp_FP_Cell_Order[i]].time = 0.0;
                
                //Right Cell
                FP_Cell[tmp_FP_Cell_Order[i]].number = 4;
                FP_Cell[tmp_FP_Cell_Order[i]].index[0] = FP_Vertex_Number-2;
                FP_Cell[tmp_FP_Cell_Order[i]].index[1] = index[1];
                FP_Cell[tmp_FP_Cell_Order[i]].index[2] = index[2];
                FP_Cell[tmp_FP_Cell_Order[i]].index[3] = FP_Vertex_Number-1;
                FP_Cell_Mig_Force[tmp_FP_Cell_Order[i]] = 0.0;
                
                if (t > Coloring_Time) FP_Cell[FP_Number].flag = 1;
                else FP_Cell[FP_Number].flag = 0;
                FP_Cell[FP_Number].divide_flag = 0;
                FP_Cell[FP_Number].time = 0.0;
                FP_Cell_Order[CellNumber_Counter++] = FP_Number; //Adding Cell Order
                
                //Calculation of Potential Energy for Edge Length
                Potential_Energy(FP_Vertex, FP_Cell, FP_Number, FP_Cell_Energy_Data);

                FP_Cell_Order[CellNumber_Counter++] = tmp_FP_Cell_Order[i];
                
                sprintf(filename, "./Data_%d/FP_Division_Data_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
                fp = fopen(filename, "a");
                fprintf(fp, "%f %d\n", t, tmp_FP_Cell_Order[i]);
                fprintf(fp, "%f %d\n", t, FP_Number);
                fclose(fp);
                
                //Adding Cell Number
                FP_Number++;
            } else {
                FP_Cell_Order[CellNumber_Counter++] = tmp_FP_Cell_Order[i];
            }
        }

        //HC Cell Division
        tmp_CellNumber = HC_Number;
        for(int i=0; i<HC_Number; i++) {
            tmp_HC_Cell_Order[i] = HC_Cell_Order[i]; }
        CellNumber_Counter = 0;
        for (int i=0; i<HC_Number; i++) {
            if ((HC_Cell_Energy_Data[tmp_HC_Cell_Order[i]] > Threshold &&
                 HC_Cell[tmp_HC_Cell_Order[i]].time > Divide_Timing)
                && HC_Number < Max_CellNumber) {

                double tmp_HC_Centroid_x = HC_CellCentroid[tmp_HC_Cell_Order[i]].x;
                double tmp_length = (tmp_HC_Centroid_x - Min_HC_Centroid_x) / (Max_HC_Centroid_x - Min_HC_Centroid_x);
 
                Vector2D Bottom_Vertex = (HC_Vertex[ HC_Cell[tmp_HC_Cell_Order[i]].index[0] ]
                                          + HC_Vertex[ HC_Cell[tmp_HC_Cell_Order[i]].index[1] ]) * 0.5;
                Vector2D Top_Vertex = (HC_Vertex[ HC_Cell[tmp_HC_Cell_Order[i]].index[2] ]
                                       + HC_Vertex[ HC_Cell[tmp_HC_Cell_Order[i]].index[3] ]) * 0.5;
                HC_Vertex[HC_Vertex_Number++] = Bottom_Vertex;
                HC_Vertex[HC_Vertex_Number++] = Top_Vertex;

                int index[4];
                for (int j=0; j<4; j++) index[j] = HC_Cell[tmp_HC_Cell_Order[i]].index[j];
                
                //Left Cell
                HC_Cell[HC_Number].number = 4;
                HC_Cell[HC_Number].index[0] = index[0];
                HC_Cell[HC_Number].index[1] = HC_Vertex_Number-2;
                HC_Cell[HC_Number].index[2] = HC_Vertex_Number-1;
                HC_Cell[HC_Number].index[3] = index[3];
                HC_Cell_Mig_Force[HC_Number] = 0.0;

                if (t > Coloring_Time) HC_Cell[tmp_HC_Cell_Order[i]].flag = 1;
                else HC_Cell[tmp_HC_Cell_Order[i]].flag = 0;
                HC_Cell[tmp_HC_Cell_Order[i]].divide_flag = 0;
                HC_Cell[tmp_HC_Cell_Order[i]].time = 0.0;
                
                //Right Cell
                HC_Cell[tmp_HC_Cell_Order[i]].number = 4;
                HC_Cell[tmp_HC_Cell_Order[i]].index[0] = HC_Vertex_Number-2;
                HC_Cell[tmp_HC_Cell_Order[i]].index[1] = index[1];
                HC_Cell[tmp_HC_Cell_Order[i]].index[2] = index[2];
                HC_Cell[tmp_HC_Cell_Order[i]].index[3] = HC_Vertex_Number-1;
                HC_Cell_Mig_Force[tmp_HC_Cell_Order[i]] = 0.0;

                if (t > Coloring_Time) HC_Cell[HC_Number].flag = 1;
                else HC_Cell[HC_Number].flag = 0;
                HC_Cell[HC_Number].divide_flag = 0;
                HC_Cell[HC_Number].time = 0.0;
                HC_Cell_Order[CellNumber_Counter++] = HC_Number; //Adding Cell Order
                
                //Calculation of Potential Energy for Edge Length
                Potential_Energy(HC_Vertex, HC_Cell, HC_Number, HC_Cell_Energy_Data);

                HC_Cell_Order[CellNumber_Counter++] = tmp_HC_Cell_Order[i];
                
                //Adding Cell Number
                HC_Number++;
            } else {
                HC_Cell_Order[CellNumber_Counter++] = tmp_HC_Cell_Order[i];
            }
        }

        Vector2D tmp_Centroid;
        for (int i=0; i<FP_Number; i++) {
            tmp_Centroid = Zero2D;
            for (int j=0; j<4; j++) { tmp_Centroid += FP_Vertex[ FP_Cell[i].index[j] ] * 0.25; }
            FP_CellCentroid[i] = tmp_Centroid;
        }
        for (int i=0; i<HC_Number; i++) {
            tmp_Centroid = Zero2D;
            for (int j=0; j<4; j++) { tmp_Centroid += HC_Vertex[ HC_Cell[i].index[j] ] * 0.25; }
            HC_CellCentroid[i] = tmp_Centroid;
        }
        
        //NOTO Cell Division
        if (i_time % Divide_Step == 0) {
            Vector2D Upper_Middle_Vertex;
            Vector2D Lower_Middle_Vertex;
            int tmp_divide_index = 1;

            Upper_Middle_Vertex = (NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[2]]
            + NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[3]]) * 0.5;
            Lower_Middle_Vertex = (NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[1]]
            + NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[0]]) * 0.5;
            NOTO_Vertex[NOTO_Vertex_Number++] = Lower_Middle_Vertex;
            NOTO_Vertex[NOTO_Vertex_Number++] = Upper_Middle_Vertex;

            NOTO_Cell[NOTO_Number].number = NOTO_Number;
            NOTO_Cell[NOTO_Number].index[0] = NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[0];
            NOTO_Cell[NOTO_Number].index[1] = NOTO_Vertex_Number-2;
            NOTO_Cell[NOTO_Number].index[2] = NOTO_Vertex_Number-1;
            NOTO_Cell[NOTO_Number].index[3] = NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[3];
            NOTO_Number++;
            
            NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[0] = NOTO_Vertex_Number-2;
            NOTO_Cell[Init_Cell_Number*2-tmp_divide_index].index[3] = NOTO_Vertex_Number-1;
        }
        
        //Calculation of Vertex Model
        //FP Cell: Time Evolution
        Nonlinear_Term(FP_Vertex, FP_Vertex_Number, FP_Cell, FP_Number, FP_Force, 1, FP_Cell_Mig_Force);
        for (int i=0; i<FP_Vertex_Number; i++) {
            Vector2D tmp_Rand_Vector = {Normal_Distribution(0, 1), Normal_Distribution(0, 1)};
            New_FP_Vertex[i] = FP_Vertex[i] + FP_Force[i] * dt
            + sigma * sqrt(dt) * tmp_Rand_Vector;
            //Movement of Left Side Vertices is controled by Control_Parameter
            if (i == FixedVertexBottom || i == FixedVertexTop) {
                New_FP_Vertex[i].x *= Control_Parameter;
            }
        }
        
        //HC Cell: Time Evolution
        Nonlinear_Term(HC_Vertex, HC_Vertex_Number, HC_Cell, HC_Number, HC_Force, 2, HC_Cell_Mig_Force);
        for (int i=0; i<HC_Vertex_Number; i++) {
            Vector2D tmp_Rand_Vector = {Normal_Distribution(0, 1), Normal_Distribution(0, 1)};
            New_HC_Vertex[i] = HC_Vertex[i] + HC_Force[i] * dt
            + sigma * sqrt(dt) * tmp_Rand_Vector;
            //Movement of Left Side Vertices is controled by Control_Parameter
            if (i == FixedVertexBottom || i == FixedVertexTop) {
                New_HC_Vertex[i].x *= Control_Parameter;
            }
        }

        //NOTO Cell: Time Evolution
        Nonlinear_Term(NOTO_Vertex, NOTO_Vertex_Number, NOTO_Cell, NOTO_Number, NOTO_Force, 0, NULL);
        for (int i=0; i<NOTO_Vertex_Number; i++) {
            Vector2D tmp_Rand_Vector = {Normal_Distribution(0, 1), Normal_Distribution(0, 1)};
            New_NOTO_Vertex[i] = NOTO_Vertex[i] + NOTO_Force[i] * dt
            + sigma * sqrt(dt) * tmp_Rand_Vector;
            //Movement of Left Side Vertices is controled by Control_Parameter
            if (i == FixedVertexBottom_NOTO || i == FixedVertexTop_NOTO) {
                New_NOTO_Vertex[i].x *= Control_Parameter;
            }
        }

        //Update of Vertices
        Vector2DArraySwap(&FP_Vertex,&New_FP_Vertex);
        Vector2DArraySwap(&HC_Vertex,&New_HC_Vertex);
        Vector2DArraySwap(&NOTO_Vertex,&New_NOTO_Vertex);

        //Cell Adhesion between NOTO Cell and Last FP-HC Cell
#ifdef Cell_Adhesion
//        FP_Vertex[FP_Cell[Init_Cell_Number-1].index[0]] = NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-1].index[3]];
//        HC_Vertex[HC_Cell[Init_Cell_Number-1].index[3]] = NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-1].index[0]];

        FP_Vertex[FP_Cell[Init_Cell_Number-1].index[1]] = NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-1].index[2]];
        HC_Vertex[HC_Cell[Init_Cell_Number-1].index[2]] = NOTO_Vertex[NOTO_Cell[Init_Cell_Number*2-1].index[1]];
#endif

        //Calculation of Potential Energy for Edge Length (Apical and Basal)
        Potential_Energy(FP_Vertex, FP_Cell, FP_Number, FP_Cell_Energy_Data);
        Potential_Energy(HC_Vertex, HC_Cell, HC_Number, HC_Cell_Energy_Data);
        
        //Migration Force: Runge-Kutta 4-4 (RK4) method
        //Variable for Runge-Kutta
        double k1, k2, k3, k4;
        double tmp_Cell_Mig_Force;
        double max_position = FP_CellCentroid[0].x;
        double min_position = FP_CellCentroid[0].x;
        for (int i=0; i<FP_Number; i++) {
            if(max_position < FP_CellCentroid[i].x) max_position = FP_CellCentroid[i].x;
            if(min_position > FP_CellCentroid[i].x) min_position = FP_CellCentroid[i].x;
        }
        for (int i=0; i<FP_Number; i++) {
            double tmp_velocity_value = Cell_Migration_Force(max_position, min_position, FP_CellCentroid[FP_Cell_Order[i]].x);
            k1 = eta * (tmp_velocity_value - FP_Cell_Mig_Force[FP_Cell_Order[i]]);
            tmp_Cell_Mig_Force = FP_Cell_Mig_Force[FP_Cell_Order[i]] + k1 * dt * 0.5;
            k2 = eta * (tmp_velocity_value - tmp_Cell_Mig_Force);
            tmp_Cell_Mig_Force = FP_Cell_Mig_Force[FP_Cell_Order[i]] + k2 * dt * 0.5;
            k3 = eta * (tmp_velocity_value - tmp_Cell_Mig_Force);
            tmp_Cell_Mig_Force = FP_Cell_Mig_Force[FP_Cell_Order[i]] + k3 * dt;
            k4 = eta * (tmp_velocity_value - tmp_Cell_Mig_Force);
            FP_Cell_Mig_Force[FP_Cell_Order[i]] = FP_Cell_Mig_Force[FP_Cell_Order[i]]
            + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
        }
        max_position = HC_CellCentroid[0].x;
        min_position = HC_CellCentroid[0].x;
        for (int i=0; i<HC_Number; i++) {
            if(max_position < HC_CellCentroid[i].x) max_position = HC_CellCentroid[i].x;
            if(min_position > HC_CellCentroid[i].x) min_position = HC_CellCentroid[i].x;
        }
        for (int i=0; i<HC_Number; i++) {
            double tmp_velocity_value = Cell_Migration_Force(max_position, min_position, HC_CellCentroid[HC_Cell_Order[i]].x);
            k1 = eta * (tmp_velocity_value - HC_Cell_Mig_Force[HC_Cell_Order[i]]);
            tmp_Cell_Mig_Force = HC_Cell_Mig_Force[HC_Cell_Order[i]] + k1 * dt * 0.5;
            k2 = eta * (tmp_velocity_value - tmp_Cell_Mig_Force);
            tmp_Cell_Mig_Force = HC_Cell_Mig_Force[HC_Cell_Order[i]] + k2 * dt * 0.5;
            k3 = eta * (tmp_velocity_value - tmp_Cell_Mig_Force);
            tmp_Cell_Mig_Force = HC_Cell_Mig_Force[HC_Cell_Order[i]] + k3 * dt;
            k4 = eta * (tmp_velocity_value - tmp_Cell_Mig_Force);
            HC_Cell_Mig_Force[HC_Cell_Order[i]] = HC_Cell_Mig_Force[HC_Cell_Order[i]]
            + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
        }
    }
    
    sprintf(filename, "mkdir Vertex_Data_%d", test_index);
    system(filename);
    
    //Vertex Data of FP Cells
    sprintf(filename, "./Vertex_Data_%d/FP_Vertex_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
    fp = fopen(filename, "w");
    fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
    for (int i=0; i<FP_Number; i++) {
        fprintf(fp, "%f\n", FP_CellCentroid[i].x);
    }
    fclose(fp);
    //Vertex Data of HC Cells
//    sprintf(filename, "./Vertex_Data_%d/HC_Vertex_%4.4f_%4.4f.dat", test_index, Velocity_Index, Migration_Parameter);
//    fp = fopen(filename, "w");
//    fprintf(fp, "%f %f\n", Velocity_Index, Migration_Parameter);
//    for (int i=0; i<HC_Number; i++) {
//        fprintf(fp, "%f\n", HC_CellCentroid[i].x);
//    }
//    fclose(fp);
    
    free(FP_Vertex);
    free(HC_Vertex);
    free(NOTO_Vertex);
    free(New_FP_Vertex);
    free(New_HC_Vertex);
    free(New_NOTO_Vertex);
    
    free(FP_Force);
    free(HC_Force);
    free(NOTO_Force);
    
    free(FP_Cell_Energy_Data);
    free(HC_Cell_Energy_Data);
    
    free(FP_Cell_Order);
    free(tmp_FP_Cell_Order);
    free(HC_Cell_Order);
    free(tmp_HC_Cell_Order);
    
    free(FP_CellCentroid);
    free(HC_CellCentroid);
    
    free(FP_Cell_Mig_Force);
    free(HC_Cell_Mig_Force);
    
    return 0;
}
