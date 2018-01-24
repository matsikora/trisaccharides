/* Written by Robert Davis to detect clashes in a structure
 * Functions modified by Oliver Grant
 * Modifications to functions are tagged with //OG 13Jun2012
 * Additions to main are not tagged
 * Additions allow multiple glycans to be built from a text file
 */
#include <algorithm>
#include <iostream>
#include <vector>
#include <mpi.h>
#include <fstream>
#include <sstream>

#include "gmml/gmml.h"

using namespace gmml;
using namespace std;
using namespace gmml::carbohydrate;

using std::list;
using std::string;

using gmml::glycam_build;
using gmml::load_parameter_file;
using gmml::load_prep_file;
using gmml::Structure;

using gmml::carbohydrate::GCBStructure;
using gmml::carbohydrate::GlycanConformationBuilder;
using gmml::carbohydrate::set_phi;


void load_files();
char check_for_clashes(const Structure& structure);
void print_rings(const vector<vector<size_t> >& rings);
CoordinateGrid<int> *init_coordinate_grid(const Structure& structure,
                                          int grid_unit);
Coordinate *get_ring_center(const Structure& structure,
                            const vector<size_t>& ring);

const double kDistanceThreshold = 2.5;

char link28='N';
int main(int argc, char *argv[]) {
    
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //enable_warnings();

    load_files();
    char clash='N';

    //string sequence = "DGalpb1-4DGlcpNAcb1-6[DGalpb1-4DGlcpNAcb1-2]DManpa1-6[DGalpb1-4[DGalpb1-4DGlcpNAcb1-6]DGlcpNAcb1-4[DGalpb1-4DGlcpNAcb1-2]DManpa1-3]DManpb1-OH";
// OG addition
   // clash=check_for_clashes(*s);
    ifstream file("GLYCANS.txt");
    string line;
    int line_number = 0;
    while (getline(file, line)) {
        if (line_number++%size != rank)
            continue;
        istringstream ss(line);
        string id;
        string glycan;
        ss >> id >> glycan;
        cout << "building glycan " << id << endl;

        Structure *s = NULL;
        try {
            s = glycam_build(glycan);
        }catch (const std::exception& exception) {
            cout << "Error on glycan " << id << ": ";
            cout << exception.what() << endl;
            continue;
        }
	
        GlycanConformationBuilder b(*s);
        ///* fast
        b.add_likely_omega_values();
        b.add_phi_value("Neu5Ac", 2, "*", 6, -60.0);
        //b.add_phi_value("Neu5Ac", 2, "*", 6, 180.0); //5% popn
        b.add_phi_value("Neu5Ac", 2, "*", 3, 180.0);
        b.add_phi_value("Neu5Ac", 2, "*", 3, -60.0);
        //Neu5Gc
        b.add_phi_value("Neu5Gc", 2, "*", 3, 180.0);
        b.add_phi_value("Neu5Gc", 2, "*", 3, -60.0);
        //b.add_phi_value("Neu5Gc", 2, "*", 6, 180.0);
        b.add_phi_value("Neu5Gc", 2, "*", 6, -60.0);
        //KDN
        //b.add_phi_value("KDN", 2, "*", 3, 180.0);
        //b.add_phi_value("KDN", 2, "*", 3, -60.0);
        //b.add_phi_value("KDN", 2, "*", 6, 180.0);
        //b.add_phi_value("KDN", 2, "*", 6, -60.0);
        //fast */ 

        int index = 0;
        list<GCBStructure*> *structures = b.build();
        list<GCBStructure*>::iterator it = structures->begin();
        bool set;
        while (it != structures->end()) {
/*******************Horrible code for 2-8 linkages******************************/
            link28='N';
            for (int i = 0; i < (*it)->residue_count(); i++) {
                string code = (*it)->residues(i)->name();
                //cout << "name=" << code << endl;
                if (code == "8SA" || code == "8SB") {
                    //STATE A
                    link28='Y';
                    set=set_phi(s,i+1,-62.0);
                    (*it)->set_dihedral(i+1, "C2", i, "O8", i, "C8", i, "H8", 1.0);
                    (*it)->set_dihedral(i, "H8", i, "C8", i, "C7", i, "H7", 81.0);
                    (*it)->set_dihedral(i, "H7", i, "C7", i, "C6", i, "H6",  -64.0);
                    //(*it)->minimize("min.in");
                    clash='N'; // reset
                    clash=check_for_clashes((**it));
                    if (clash=='N') {
                        //cout << "building link28 torsion " << (*it)->name() << endl;
                        (*it)->print_amber_top_file("MD_Files/" + id + "-" + "A-" + (*it)->name() + ".prmtop");
                        (*it)->print_coordinate_file("MD_Files/" + id + "-" + "A-" + (*it)->name() + ".rst7");
                        (*it)->print_pdb_file("pdb_Files/" + id + "-" + "A-" + (*it)->name() + ".pdb");
                    }
                    //STATE B
                    set=set_phi(s,i+1,-76.0);
                    (*it)->set_dihedral(i+1, "C2", i, "O8", i, "C8", i, "H8", 16.0);
                    (*it)->set_dihedral(i, "H8", i, "C8", i, "C7", i, "H7", -58.0);
                    (*it)->set_dihedral(i, "H7", i, "C7", i, "C6", i, "H6",  -60.0);
                    //(*it)->minimize("min.in");
                    clash='N'; // reset
                    clash=check_for_clashes((**it));
                    if (clash=='N') {
                        //cout << "building link28 torsion " << (*it)->name() << endl;
                        (*it)->print_amber_top_file("MD_Files/" + id + "-" + "B-" + (*it)->name() + ".prmtop");
                        (*it)->print_coordinate_file("MD_Files/" + id + "-" + "B-" + (*it)->name() + ".rst7");
                        (*it)->print_pdb_file("pdb_Files/" + id + "-" + "B-" + (*it)->name() + ".pdb");
                    }
                    //STATE C
                    set=set_phi(s,i+1,-76.0);
                    (*it)->set_dihedral(i+1, "C2", i, "O8", i, "C8", i, "H8", -53.0);
                    (*it)->set_dihedral(i, "H8", i, "C8", i, "C7", i, "H7", 74.0);
                    (*it)->set_dihedral(i, "H7", i, "C7", i, "C6", i, "H6",  -64.0);
                    //(*it)->minimize("min.in");
                    clash='N'; // reset
                    clash=check_for_clashes((**it));
                    if (clash=='N') {
                        //cout << "building link28 torsion " << (*it)->name() << endl;
                        (*it)->print_amber_top_file("MD_Files/" + id + "-" + "C-" + (*it)->name() + ".prmtop");
                        (*it)->print_coordinate_file("MD_Files/" + id + "-" + "C-" + (*it)->name() + ".rst7");
                        (*it)->print_pdb_file("pdb_Files/" + id + "-" + "C-" + (*it)->name() + ".pdb");
                    }
                    //STATE D
                    set=set_phi(s,i+1,-50.0);
                    (*it)->set_dihedral(i+1, "C2", i, "O8", i, "C8", i, "H8", 63.0);
                    (*it)->set_dihedral(i, "H8", i, "C8", i, "C7", i, "H7", -178.0);
                    (*it)->set_dihedral(i, "H7", i, "C7", i, "C6", i, "H6",  -60.0);
                    //(*it)->minimize("min.in");
                    clash='N'; // reset
                    clash=check_for_clashes((**it));
                    if (clash=='N') {
                        //cout << "building link28 torsion " << (*it)->name() << endl;
                        (*it)->print_amber_top_file("MD_Files/" + id + "-" + "D-" + (*it)->name() + ".prmtop");
                        (*it)->print_coordinate_file("MD_Files/" + id + "-" + "D-" + (*it)->name() + ".rst7");
                        (*it)->print_pdb_file("pdb_Files/" + id + "-" + "D-" + (*it)->name() + ".pdb");
                    }
                    //STATE E
                    set=set_phi(s,i+1,-164.0);
                    (*it)->set_dihedral(i+1, "C2", i, "O8", i, "C8", i, "H8", -15.0);
                    (*it)->set_dihedral(i, "H8", i, "C8", i, "C7", i, "H7", 77.0);
                    (*it)->set_dihedral(i, "H7", i, "C7", i, "C6", i, "H6",  -57.0);
                    //(*it)->minimize("min.in");
                    clash='N'; // reset
                    clash=check_for_clashes((**it));
                    if (clash=='N') {
                        //cout << "building link28 torsion " << (*it)->name() << endl;
                        (*it)->print_amber_top_file("MD_Files/" + id + "-" + "E-" + (*it)->name() + ".prmtop");
                        (*it)->print_coordinate_file("MD_Files/" + id + "-" + "E-" + (*it)->name() + ".rst7");
                        (*it)->print_pdb_file("pdb_Files/" + id + "-" + "E-" + (*it)->name() + ".pdb");
                    }
                    //STATE F
                    set=set_phi(s,i+1,-52.0);
                    (*it)->set_dihedral(i+1, "C2", i, "O8", i, "C8", i, "H8", 49.0);
                    (*it)->set_dihedral(i, "H8", i, "C8", i, "C7", i, "H7", 81.0);
                    (*it)->set_dihedral(i, "H7", i, "C7", i, "C6", i, "H6",  -61.0);
                    //(*it)->minimize("min.in");
                    clash='N'; // reset
                    clash=check_for_clashes((**it));
                    if (clash=='N') {
                        //cout << "building link28 torsion " << (*it)->name() << endl;
                        (*it)->print_amber_top_file("MD_Files/" + id + "-" + "F-" + (*it)->name() + ".prmtop");
                        (*it)->print_coordinate_file("MD_Files/" + id + "-" + "F-" + (*it)->name() + ".rst7");
                        (*it)->print_pdb_file("pdb_Files/" + id + "-" + "F-" + (*it)->name() + ".pdb");
                    }
                }
            }
/*******************END horrible code for 2-8 linkages***************************/
            if (link28=='N') {
                //(*it)->minimize("min.in");
                clash='N'; // reset 
                clash=check_for_clashes((**it));
                if (clash=='N') {
               // cout << "building nonlink28 torsion " << (*it)->name() << endl;
                    (*it)->print_amber_top_file("MD_Files/" + id + "-" + (*it)->name() + ".prmtop");
                    (*it)->print_coordinate_file("MD_Files/" + id + "-" + (*it)->name() + ".rst7");
		    (*it)->print_pdb_file("pdb_Files/" + id + "-" + (*it)->name() + ".pdb");
                }
            }
            ++it;
            ++index;
        }
    }
    MPI_Finalize();
}

// You need to modify this.
void load_files() {
    add_path("../../");
    load_parameter_file("param_files/parm99.dat.mod");
    load_parameter_file("param_files/Glycam_06g.dat");
    load_prep_file("prep_files/Glycam_06.prep");
    load_prep_file("prep_files/Neu5Gc_a_06.prep");
    load_prep_file("prep_files/sulfate.prep");
    load_prep_file("prep_files/Neu5Ac_b_06.prep");
}

char check_for_clashes(const Structure& structure) {
    char clash='N'; //OG 13Jun2012
    // This returns a list of rings, each ring being a list of atom indices.
    // There should probably be a function in Structure called get_rings() that
    // makes the call to Graph::get_cycles(), but for now you can do this.
    // (size_t is just an unsigned integer type.)
    vector<vector<size_t> > *rings = structure.bonds()->get_cycles();
    //print_rings(*rings);

    // You can read about this in coordinate_grid.h. It's not the fanciest way
    // to do spatial indexing, but it scales well.
    CoordinateGrid<int> *grid = init_coordinate_grid(structure,
                                                     kDistanceThreshold);

    // Vertices are residue indices, there's an edge between two vertices if
    // the corresponding residues are linked.
    Graph *link_graph = structure.get_link_graph();

    for (int i = 0; i < rings->size(); i++) {
        Coordinate *center = get_ring_center(structure, (*rings)[i]);
        vector<int> *adjacent_cells = grid->retrieve_adjacent_cells(*center);

        // Since the grid stores stuff in cubes and we only want coordinates
        // within a sphere, we need to filter the atoms in the adjacent cells.
        vector<int> close_atoms;
        for (int j = 0; j < adjacent_cells->size(); j++) {
            int atom_index = (*adjacent_cells)[j];
            const Atom *close_atom = structure.atoms(atom_index);
            if (measure(close_atom->coordinate(), *center) <=
                    kDistanceThreshold)
                close_atoms.push_back(atom_index);
        }
        delete adjacent_cells;

        // We'll say the residue of the first atom in the ring is the residue
        // of the ring.
        int this_residue = structure.get_residue_index((*rings)[i][0]);

        const Graph::AdjList& adjacent_residues =
                link_graph->edges(this_residue);

        for (int j = 0; j < close_atoms.size(); j++) {
            int that_residue = structure.get_residue_index(close_atoms[j]);
            if (that_residue == this_residue)
                continue;
            // If that_residue is a residue that's not adjacent to this_residue,
            // we've got a clash.
            if (find(adjacent_residues.begin(), adjacent_residues.end(),
                     that_residue) == adjacent_residues.end()) {
                const Atom *close_atom = structure.atoms(close_atoms[j]);
                double distance = measure(*center, close_atom->coordinate());
		clash='Y'; //OG 13Jun2012
                j=close_atoms.size(); // Finish looping if find even one clash //OG 13Jun2012
                i=rings->size(); // Finish looping if find even one clash //OG 13Jun2012
                /*cout << "atom " << close_atom->name() << " in residue " <<
                       that_residue << " is too close to the ring center " <<
                        "in residue " << this_residue <<
                        " (distance: " << distance << ")" << endl; 
                */
                cout << "Clash detected " << endl;
                 
            }
        }
        delete center;
    }
    delete link_graph;
    delete grid;
    delete rings;
    return clash; //OG 13Jun2012
}

void print_rings(const vector<vector<size_t> >& rings) {
    cout << "Rings: " << endl;
    for (int i = 0; i < rings.size(); i++) {
        for (int j = 0; j < rings[i].size(); j++) {
            cout << rings[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

CoordinateGrid<int> *init_coordinate_grid(const Structure& structure,
                                          int grid_unit) {
    CoordinateGrid<int> *grid = new CoordinateGrid<int>(grid_unit);
    vector<size_t> *residue_index_table =
            structure.get_residue_index_table();
    for (int i = 0; i < structure.size(); i++) {
        grid->insert(structure.atoms(i)->coordinate(), i);
    }
    delete residue_index_table;
    return grid;
}

// This just averages the ring coordinates. I don't know how good of a center
// this'll give you.
Coordinate *get_ring_center(const Structure& structure,
                            const vector<size_t>& ring) {
    double x_sum = 0.0;
    double y_sum = 0.0;
    double z_sum = 0.0;

    int ring_size = ring.size();

    for (int i = 0; i < ring_size; i++) {
        const Coordinate& c = structure.atoms(ring[i])->coordinate();
        x_sum += c.x;
        y_sum += c.y;
        z_sum += c.z;
    }
    
    return new Coordinate(x_sum/ring_size, y_sum/ring_size, z_sum/ring_size);
}
