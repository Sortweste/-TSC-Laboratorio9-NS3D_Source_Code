#include <fstream>
#include "string.h"

void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    if(nlines==DOUBLELINE) file >> line;

    for(int i=0;i<n;i++){
        switch(mode){
            case INT_FLOAT:
                int e0; float condition_value;
                file >> e0 >> condition_value;
                item_list[i].setValues(0,0,0,0,0,e0,condition_value,0);
                break;
            case INT_FLOAT_FLOAT_FLOAT:
                int e; float coord1, coord2, coord3;
                file >> e >> coord1 >> coord2 >> coord3;
                item_list[i].setValues(e,coord1,coord2,coord3,0,0,0,0);
                break;
            case INT_INT_INT_INT_INT:
                int element,node1,node2,node3,node4;
                file >> element >> node1 >> node2 >> node3 >> node4;
                item_list[i].setValues(element,0,0,0,node1,node2,node3,node4);
                break;
            }
    }
}

void correctConditions(int n,condition *list,int *indices){
    for(int i=0;i<n;i++)
        indices[i] = list[i].getNode1();

    for(int i=0;i<n-1;i++){
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}

void addExtension(char *newfilename,char *filename, string extension){
    int ori_length = strlen(filename);
    int ext_length = extension.length();
    int i;
    for(i=0;i<ori_length;i++)
        newfilename[i] = filename[i];
    for(i=0;i<ext_length;i++)
        newfilename[ori_length+i] = extension[i];
    newfilename[ori_length+i] = '\0';
}

void correctIndices(int n,condition* list,int total){
    for(int i=0;i<n;i++)
        list[i].setNode1(list[i].getNode1()+total);
}

void addArray(int* index, int n,condition* list, condition* listf){
    for(int i=0;i<n;i++){
        listf[*index] = list[i];
        (*index)++;
    }
}

void fusionDirichlet(int n1,condition* list1,int n2,condition* list2,int n3,condition* list3,int n4,condition* list4,condition* list){
    int index = 0;
    addArray(&index,n1,list1,list);
    addArray(&index,n2,list2,list);
    addArray(&index,n3,list3,list);
    addArray(&index,n4,list4,list);
}

void leerMallaYCondiciones(mesh &m,char *filename){
    char inputfilename[150];
    ifstream file;
    float u_bar,nu,rho,f_x,f_y,f_z;
    int nnodes, neltos, ndirich_u, ndirich_v, ndirich_w, ndirich_p;
    condition *dirichlet_u, *dirichlet_v, *dirichlet_w, *dirichlet_p;

    addExtension(inputfilename,filename,".dat");
    file.open(inputfilename);

    file >> u_bar >> nu >> rho >> f_x >> f_y >> f_z;
    
    file >> nnodes >> neltos >> ndirich_u >> ndirich_v >> ndirich_w >> ndirich_p;

    m.setParameters(u_bar,nu,rho,f_x,f_y,f_z);
    m.setSizes(nnodes,neltos,ndirich_u+ndirich_v+ndirich_w+ndirich_p);
    m.createData();

    dirichlet_u = new condition[ndirich_u];
    dirichlet_v = new condition[ndirich_v];
    dirichlet_w = new condition[ndirich_w];
    dirichlet_p = new condition[ndirich_p];

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT_FLOAT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT_INT_INT,m.getElements());
    obtenerDatos(file,DOUBLELINE,ndirich_u,INT_FLOAT,dirichlet_u);
    obtenerDatos(file,DOUBLELINE,ndirich_v,INT_FLOAT,dirichlet_v);
    obtenerDatos(file,DOUBLELINE,ndirich_w,INT_FLOAT,dirichlet_w);
    obtenerDatos(file,DOUBLELINE,ndirich_p,INT_FLOAT,dirichlet_p);

    file.close();

    correctIndices(ndirich_v,dirichlet_v,nnodes);
    correctIndices(ndirich_w,dirichlet_w,2*nnodes);
    correctIndices(ndirich_p,dirichlet_p,3*nnodes);

    fusionDirichlet(ndirich_u,dirichlet_u,ndirich_v,dirichlet_v,ndirich_w,dirichlet_w,ndirich_p,dirichlet_p,m.getDirichlet());
    
    correctConditions(ndirich_u+ndirich_v+ndirich_w+ndirich_p,m.getDirichlet(),m.getDirichletIndices());
}

bool findIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return true;
    return false;
}

int getIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return i;
    return -1;
}

int getIndex(int v, int s, Vector vec){
    for(int i=0;i<s;i++)
        if(vec.at(i)==v) return i;
    return -1;
}

int *createNonDirichletIndices(int nn,int nd,int *dirich_indices){
    int *ndi = new int[4*nn-nd];
    int pos = 0;
    for(int i=1;i<=4*nn;i++)
        if(!findIndex(i,nd,dirich_indices)){
            ndi[pos] = i;
            pos++;
        }
    return ndi;
}

void writeResults(mesh m,Vector T,char *filename){
    char outputfilename[150];

    int nn = m.getSize(NODES);
    int nd = m.getSize(DIRICHLET);
    int *dirich_indices = m.getDirichletIndices();
    int *non_dirich_indices = createNonDirichletIndices(nn,nd,dirich_indices);
    condition *dirich = m.getDirichlet();
    ofstream file;

    addExtension(outputfilename,filename,".post.res");
    file.open(outputfilename);

    file << "GiD Post Results File 1.0\n";

    file << "Result \"Velocity\" \"Load Case 1\" 1 Vector OnNodes\nComponentNames \"u\" \"v\" \"w\"\nValues\n";

    for(int i=0;i<nn;i++){
        int d_index = getIndex(i+1,nd,dirich_indices);
        if(d_index != -1)
            file << i+1 << " " << dirich[d_index].getValue() << " ";
        else{
            int T_index = getIndex(i+1,4*nn-nd,non_dirich_indices);
            file << i+1 << " " << T.at(T_index) << " ";
        }
        d_index = getIndex(i+1+nn,nd,dirich_indices);
        if(d_index != -1)
            file << dirich[d_index].getValue() << "\n";
        else{
            int T_index = getIndex(i+1+nn,4*nn-nd,non_dirich_indices);
            file << T.at(T_index) << "\n";
        }
    }

    file << "End values\n";

    file << "\nResult \"Pressure\" \"Load Case 1\" 1 Scalar OnNodes\nComponentNames \"p\"\nValues\n";

    for(int i=0;i<nn;i++){
        int d_index = getIndex(i+1+2*nn,nd,dirich_indices);
        if(d_index != -1)
            file << i+1 << " " << dirich[d_index].getValue() << "\n";
        else{
            int T_index = getIndex(i+1+2*nn,4*nn-nd,non_dirich_indices);
            file << i+1 << " " << T.at(T_index) << "\n";
        }
    }

    file << "End values\n";

    file.close();
}
