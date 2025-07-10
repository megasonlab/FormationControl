//Function.h
#ifndef _FUNCTION_H
#define _FUNCTION_H

//Structure of Cell
typedef struct{
    int number;
    int index[4];
    int flag;
    int divide_flag;
    double time;
} Cell;

double Uniform_Random_Number(){
    return genrand_real3();
}
double Normal_Distribution(double mu, double sigma){
    return mu + sigma * sqrt(-2.0 * log(Uniform_Random_Number())) * sin(2.0 * Pi * Uniform_Random_Number());
}

Vector2D NormalizedVector(Vector2D a) {
    return a / ~a;
}

double Triangle_Area(Vector2D a, Vector2D b, Vector2D g) {
    return sqrt(~(a-g)*~(a-g)*~(b-g)*~(b-g) - ((a-g)*(b-g))*((a-g)*(b-g))) * 0.5;
}

Vector2D Function_Area(Vector2D a, Vector2D b, Vector2D c) {
    Vector3D tmp_a = {a.x,a.y,0.0};
    Vector3D tmp_b = {b.x,b.y,0.0};
    Vector3D tmp_c = {c.x,c.y,0.0};
    Vector3D tmp_ans = ((tmp_a^tmp_b)^tmp_c) / (2.0 * ~(tmp_a^tmp_b));
    Vector2D ans = {tmp_ans.x, tmp_ans.y};
    return ans;
}

void Potential_Energy(Vector2D *Vertex, Cell *Cell_Data, int CellNumber, double *Energy) {    
    for (int i=0; i<CellNumber; i++) {
        double tmp_Basal_Length = ~(Vertex[ Cell_Data[i].index[0] ] - Vertex[ Cell_Data[i].index[1] ]);
        double tmp_Apical_Length = ~(Vertex[ Cell_Data[i].index[2] ] - Vertex[ Cell_Data[i].index[3] ]);
        Energy[i] = b * ( tmp_Basal_Length * tmp_Basal_Length + tmp_Apical_Length * tmp_Apical_Length ) * 0.5;
    }
}

void Nonlinear_Term(Vector2D *Vertex, int VertexNumber, Cell *Cell_Data, int CellNumber,
                    Vector2D *Vertex_Force, int CellType, double *Cell_Mig_Force) {
    
    double Ideal_Area = 0.0;
    if (CellType == 0) Ideal_Area = Ideal_NOTO_Area;
    if (CellType == 1) Ideal_Area = Ideal_FP_Area;
    if (CellType == 2) Ideal_Area = Ideal_HC_Area;

    //Initialized Variable (Vertex_Force)
    for (int i=0; i<VertexNumber; i++) { Vertex_Force[i] = {0.0, 0.0}; }

    for (int i=0; i<CellNumber; i++) {
        //Cell Volume Calculation
        Vector2D tmp_Centroid = Zero2D;
        for (int j=0; j<4; j++) { tmp_Centroid += Vertex[ Cell_Data[i].index[j] ] * 0.25; }
        double tmp_Area = 0.0;
        for (int j=0; j<4; j++) {
            Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
            Vector2D Next_Vertex = Vertex[ Cell_Data[i].index[(j+1)%4] ];
            tmp_Area += Triangle_Area(Next_Vertex, Current_Vertex, tmp_Centroid);
        }

        
        for (int j=0; j<4; j++) {
            Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
            Vector2D Prev_Vertex = Vertex[ Cell_Data[i].index[(j+3)%4] ];
            Vector2D Next_Vertex = Vertex[ Cell_Data[i].index[(j+1)%4] ];
            
            //Cell Volume Conservation
            Vertex_Force[ Cell_Data[i].index[j] ]
            += a * Function_Area(Next_Vertex-Current_Vertex, Prev_Vertex-Current_Vertex, Next_Vertex-Prev_Vertex)
            * (Ideal_Area - tmp_Area) * (-1.0);

            //Cell Boundary Contraction
            Vertex_Force[ Cell_Data[i].index[j] ]
            += b * NormalizedVector(Current_Vertex - Next_Vertex)
            * ~(Current_Vertex - Next_Vertex) * (-1.0);
            Vertex_Force[ Cell_Data[i].index[j] ]
            += b * NormalizedVector(Current_Vertex - Prev_Vertex)
            * ~(Current_Vertex - Prev_Vertex) * (-1.0);
        }
        
        //Adhesion of Basal Edge
        //NOTO Cells
        double Membrane_Length = 0.0;
        if (CellType == 0) {
            Membrane_Length = NOTO_Height * 0.5;
            for (int j=0; j<=1; j++) {
                Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
                Vector2D Membrane_Vertex = {Current_Vertex.x, 0.0};
                Vector2D tmp_Vector = Membrane_Vertex - Current_Vertex;
                double tmp_Length = ~tmp_Vector;
                if (1.0e-15 < tmp_Length && tmp_Length < (2.0*Membrane_Length + 0.1) * 0.5) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
                if ((2.0*Membrane_Length + 0.1) * 0.5 < tmp_Length && tmp_Length < Membrane_Length + 0.1) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
            }
        }
        if (CellType == 0) {
            Membrane_Length = NOTO_Height * 0.5;
            for (int j=2; j<=3; j++) {
                Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
                Vector2D Membrane_Vertex = {Current_Vertex.x, 0.0};
                Vector2D tmp_Vector = Membrane_Vertex - Current_Vertex;
                double tmp_Length = ~tmp_Vector;
                if (1.0e-15 < tmp_Length && tmp_Length < (2.0*Membrane_Length + 0.1) * 0.5) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
                if ((2.0*Membrane_Length + 0.1) * 0.5 <= tmp_Length && tmp_Length < Membrane_Length + 0.1) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
            }
        }
        //FP Cells
        if (CellType == 1) {
            Membrane_Length = NOTO_Height * 0.5;
            for (int j=0; j<=1; j++) {
                Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
                Vector2D Membrane_Vertex = {Current_Vertex.x, 0.0};
                Vector2D tmp_Vector = Membrane_Vertex - Current_Vertex;
                double tmp_Length = ~tmp_Vector;
                if (1.0e-15 < tmp_Length && tmp_Length < (2.0*Membrane_Length + 0.1) * 0.5) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
                if ((2.0*Membrane_Length + 0.1) * 0.5 <= tmp_Length && tmp_Length < Membrane_Length + 0.1) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
            }
        }
        if (CellType == 1) {
            Membrane_Length = NOTO_Height * 0.5 + FP_Height;
            for (int j=2; j<=3; j++) {
                Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
                Vector2D Membrane_Vertex = {Current_Vertex.x, 0.0};
                Vector2D tmp_Vector = Membrane_Vertex - Current_Vertex;
                double tmp_Length = ~tmp_Vector;
                if (tmp_Length > Membrane_Length)
                Vertex_Force[ Cell_Data[i].index[j] ]
                    += Apical_Repulsion_Parameter * NormalizedVector(Membrane_Vertex - Current_Vertex) 
                    * (tmp_Length - Membrane_Length) * (1.0);
            }
        }
        //HC Cells
        if (CellType == 2) {
            Membrane_Length = NOTO_Height * 0.5;
            for (int j=2; j<=3; j++) {
                Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
                Vector2D Membrane_Vertex = {Current_Vertex.x, 0.0};
                Vector2D tmp_Vector = Membrane_Vertex - Current_Vertex;
                double tmp_Length = ~tmp_Vector;
                if (1.0e-15 < tmp_Length && tmp_Length < (2.0*Membrane_Length + 0.1) * 0.5) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
                if ((2.0*Membrane_Length + 0.1) * 0.5 <= tmp_Length && tmp_Length < Membrane_Length + 0.1) {
                    Vertex_Force[ Cell_Data[i].index[j] ]
                    += Basal_Adhesion_Parameter * NormalizedVector(tmp_Vector) * (tmp_Length - Membrane_Length) * (1.0);
                }
            }
        }
        if (CellType == 2) {
            Membrane_Length = NOTO_Height * 0.5 + FP_Height;
            for (int j=0; j<=1; j++) {
                Vector2D Current_Vertex = Vertex[ Cell_Data[i].index[j] ];
                Vector2D Membrane_Vertex = {Current_Vertex.x, 0.0};
                Vector2D tmp_Vector = Membrane_Vertex - Current_Vertex;
                double tmp_Length = ~tmp_Vector;
                if (tmp_Length > Membrane_Length)
                Vertex_Force[ Cell_Data[i].index[j] ]
                    += Apical_Repulsion_Parameter * NormalizedVector(Membrane_Vertex - Current_Vertex)
                    * (tmp_Length - Membrane_Length) * (1.0);
            }
        }
    }
    
    //External Force for Cell Migration
    if (CellType == 1) {
        for (int i=0; i<=CellNumber-1; i++) {
            Vertex_Force[ Cell_Data[ i ].index[1] ]
            += {Cell_Mig_Force[ i ], 0.0};
        }
    }
   if (CellType == 2) {
       for (int i=0; i<=CellNumber-1; i++) {
           Vertex_Force[ Cell_Data[ i ].index[2] ]
           += {Cell_Mig_Force[ i ], 0.0};
       }
   }
}

double Cell_Migration_Force(double max, double min, double x) {
   return Migration_Parameter * pow((x - min)/(max - min), Velocity_Index);
}

void DoubleArraySwap(double **u, double **oldu) {
    double *tmp;
    tmp = *u;
    *u = *oldu;
    *oldu = tmp;
}

void Vector2DArraySwap(Vector2D **u, Vector2D **oldu) {
    Vector2D *tmp;
    tmp = *u;
    *u = *oldu;
    *oldu = tmp;
}

#endif
