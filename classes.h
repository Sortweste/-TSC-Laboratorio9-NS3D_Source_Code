enum lines {NOLINE,SINGLELINE,DOUBLELINE};
enum modes {NOMODE,INT_FLOAT,INT_FLOAT_FLOAT_FLOAT,INT_INT_INT_INT_INT};
enum parameters {ADJECTIVE_VELOCITY,DYNAMIC_VISCOSITY,DENSITY,EXTERNAL_FORCE_X,EXTERNAL_FORCE_Y, EXTERNAL_FORCE_Z};
enum sizes {NODES,ELEMENTS,DIRICHLET};
enum coords {EQUIS,YE,ZETA};

class item{
    protected:
        int id;
        float x;
        float y;
        float z;
        int node1;
        int node2;
        int node3;
        int node4;
        float value;
    public:
        void setId(int identifier) {
            id = identifier;
        }

        void setX(float x_coord) {
            x = x_coord;
        }

        void setY(float y_coord) {
            y = y_coord;
        }

        void setZ(float z_coord) {
            z = z_coord;
        }

        void setNode1(int node_1) {
            node1 = node_1;
        }

        void setNode2(int node_2) {
            node2 = node_2;
        }

        void setNode3(int node_3) {
            node3 = node_3;
        }

        void setNode4(int node_4) {
            node4 = node_4;
        }

        void setValue(float value_to_assign) {
            value = value_to_assign;
        }

        int getId() {
            return id;
        }

        float getX() {
            return x;
        }

        float getY() {
            return y;
        }

        float getZ() {
            return z;
        }

        int getNode1() {
            return node1;
        }

        int getNode2() {
            return node2;
        }

        int getNode3() {
            return node3;
        }

        int getNode4() {
            return node4;
        }

        float getValue() {
            return value;
        }

        virtual void setValues(int a,float b,float c,float d,int e,int f,float g, int h)=0;


};

class node: public item{

    public:
        void setValues(int a,float b,float c,float d,int e,int f,float g, int h){
            id = a;
            x = b;
            y = c;
            z = d;
        }

};

class element: public item{

    public:
        void setValues(int a,float b,float c,float d,int e,int f,float g, int h){
            id = a;
            node1 = e;
            node2 = f;
            node3 = g;
            node4 = h;
        }

};

class condition: public item{

    public:

        void setValues(int a,float b,float c,float d,int e,int f,float g, int h){
            node1 = f;
            value = g;
        }

};

class mesh{
        float parameters[6];
        int sizes[3];
        node *node_list;
        element *element_list;
        int *indices_dirich;
        condition *dirichlet_list;        
    public:
        void setParameters(float u_bar,float nu, float rho, float f_x, float f_y, float f_z){
            parameters[ADJECTIVE_VELOCITY]=u_bar;
            parameters[DYNAMIC_VISCOSITY]=nu;
            parameters[DENSITY]=rho;
            parameters[EXTERNAL_FORCE_X]=f_x;
            parameters[EXTERNAL_FORCE_Y]=f_y;
            parameters[EXTERNAL_FORCE_Z]=f_z;
        }
        void setSizes(int nnodes,int neltos,int ndirich){
            sizes[NODES] = nnodes;
            sizes[ELEMENTS] = neltos;
            sizes[DIRICHLET] = ndirich;
        }
        int getSize(int s){
            return sizes[s];
        }
        float getParameter(int p){
            return parameters[p];
        }
        void createData(){
            node_list = new node[sizes[NODES]];
            element_list = new element[sizes[ELEMENTS]];
            indices_dirich = new int[sizes[DIRICHLET]];
            dirichlet_list = new condition[sizes[DIRICHLET]];
        }
        node* getNodes(){
            return node_list;
        }
        element* getElements(){
            return element_list;
        }
        int* getDirichletIndices(){
            return indices_dirich;
        }
        condition* getDirichlet(){
            return dirichlet_list;
        }
        node getNode(int i){
            return node_list[i];
        }
        element getElement(int i){
            return element_list[i];
        }
        condition getCondition(int i, int type){
            return dirichlet_list[i];
        }
};
