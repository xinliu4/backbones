#include "H5Cpp.h"
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <ctime>
#include <map>
#include <set>
#include <utility>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "GenericIO.h"
#include "GenericIODefinitions.hpp"
#include "GenericIOMPIReader.h"
#include "GenericIOPosixReader.h"
#include "GenericIOReader.h"
#include "GenericIOUtilities.h"

using namespace std;
using namespace gio;
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

int step_ind;

// Check circular velocities with https://cds.cern.ch/record/285940/files/9508025.pdf, table 1.
// For halos ~ 2e13 (60 000 DeltaQ particles), V200 should be ~360 km/s.
// For halos ~ 3e11 (900 DeltaQ particles), V200 should be ~90 km/s.

int subfiles_count = 8;
int total_hdf5_files = 256;

int min_sod_count = 0;
int min_z0_count_threshold = 2000; // Only merger trees whose root fof mass are at least this massive will be kept; 500 to have an sod profile built at z=0, 2000 to have a concentration except fof != sod
int numbins = 20;
const float G_gravity = 4.3e-6; // Units of kpc km^2 s^-2 Msun^-1
//const float particlemass = 3.4023e+8; // Delta Quadrant mass
const float particlemass = 1.14828e+9; // Alpha Quadrant mass

int catalog_steps[100] = {42, 43, 44, 45, 46, 48, 49, 50, 52, 53, 54, 56, 57, 59, 60, 62, 63, 65, 67, 68, 70, 72, 74, 76, 77, 79, 81, 84, 86, 88, 90, 92, 95, 97, 100, 102, 105, 107, 110, 113, 116, 119, 121, 124, 127, 131, 134, 141, 144, 148, 151, 155, 159, 163, 167, 171, 176, 180, 184, 189, 194, 198, 203, 208, 213, 219, 224, 230, 235, 241, 247, 253, 259, 266, 272, 279, 286, 293, 300, 307, 315, 323, 331, 338, 347, 355, 365, 373, 382, 392, 401, 411, 421, 432, 442, 453, 464, 475, 487, 499};
const int count_steps = 100;
const int final_timestep = 499;

// hdf5 group names
const H5std_string nodeIndex_name("/forestHalos/nodeIndex");
const H5std_string descendentIndex_name("/forestHalos/descendentIndex");
const H5std_string haloTag_name("/forestHalos/haloTag");
const H5std_string nodeMass_name("/forestHalos/nodeMass");
const H5std_string timestep_name("/forestHalos/timestep");

const H5std_string firstNode_name("/forestIndex/firstNode");
const H5std_string numberOfNodes_name("/forestIndex/numberOfNodes");

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

int main(int argc, char** argv)
{

    // -------------------------------------------------------
    // ----------------- prep inputs/outputs -----------------

    int hdf5_file_ind = atoi(argv[1]);

    MPI_Init( &argc, &argv );
    int Rank, NRanks;
    MPI_Comm_rank( MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &NRanks);

    // pass out hdf5 files to ranks
    int hdf5_file_edges[subfiles_count+1];
    for(int i=0; i<subfiles_count+1; i++)
    {
        hdf5_file_edges[i] = (total_hdf5_files/subfiles_count)*i;
        cout << hdf5_file_edges[i] << endl;
    }

    int hdf5_file_0 = hdf5_file_edges[hdf5_file_ind];
    int hdf5_file_final = hdf5_file_edges[hdf5_file_ind+1];

    cout << hdf5_file_0 << " " << hdf5_file_final << endl;

    clock_t begin = clock();

    // decalre output files
    stringstream tmp_string_outfile;
    tmp_string_outfile << "AlphaQ." << hdf5_file_0 << "_" << hdf5_file_final << ".haloTags_c_vpeak_vr200";
    ofstream lightcone_c(tmp_string_outfile.str().c_str());
    lightcone_c << "root ID\tstep\thalo ID\tM200 (counts)\tcfit\tvpeak (km/s)\tv200 (km/s)\tM(r_s)\n";

    vector<int64_t> root_index_vec;
    vector<int> step_vec;
    vector<int64_t> hid_vec;
    map< pair<int, int64_t>, bool > backbone_halos_map;
    pair<int, int64_t> keykey_backbone_halos;
    stringstream tmp_string;
    

    // -------------------------------------------------------
    // ----------------- build backbones ---------------------

    // loop over each hdf5 file distributed to the current rank
    for(int hdf5_file_no=hdf5_file_0; hdf5_file_no<hdf5_file_final; hdf5_file_no++)
    {
        cout << hdf5_file_no << endl;
        cout << "elapsed time " << double(clock()-begin) / CLOCKS_PER_SEC << " s" << endl;

        // read next file
        tmp_string << "/projects/DarkUniverse_esp/childh/AlphaQ_mergertrees/AlphaQ." << hdf5_file_no << ".hdf5";
        H5std_string file_name(tmp_string.str().c_str());
        H5File f( file_name, H5F_ACC_RDONLY );
        tmp_string.str("");

        // Open datasets
        DataSet nodeIndex_set = f.openDataSet(nodeIndex_name);
        DataSet descendentIndex_set = f.openDataSet(descendentIndex_name);
        DataSet haloTags_set = f.openDataSet(haloTag_name);
        DataSet nodeMasses_set = f.openDataSet(nodeMass_name);
        DataSet timesteps_set = f.openDataSet(timestep_name);
        DataSet firstNodes_set = f.openDataSet(firstNode_name);
        DataSet numberOfNodes_set = f.openDataSet(numberOfNodes_name);

        DataSpace nodeIndex_space = nodeIndex_set.getSpace();
        DataSpace descendentIndex_space = descendentIndex_set.getSpace();
        DataSpace haloTags_space = haloTags_set.getSpace();
        DataSpace nodeMasses_space = nodeMasses_set.getSpace();
        DataSpace timesteps_space = timesteps_set.getSpace();
        DataSpace firstNodes_space = firstNodes_set.getSpace();
        DataSpace numberOfNodes_space = numberOfNodes_set.getSpace();

        hsize_t rank;
        hsize_t dims[2];

        rank = nodeIndex_space.getSimpleExtentDims(dims, NULL);
        int64_t *nodeIndex;
        nodeIndex = (int64_t *)malloc(sizeof(int64_t)*dims[0]);
        nodeIndex_set.read(nodeIndex, PredType::NATIVE_INT64, nodeIndex_space);

        rank = descendentIndex_space.getSimpleExtentDims(dims, NULL);
        int64_t *descendentIndex;
        descendentIndex = (int64_t *)malloc(sizeof(int64_t)*dims[0]);
        descendentIndex_set.read(descendentIndex, PredType::NATIVE_INT64, descendentIndex_space);

        rank = haloTags_space.getSimpleExtentDims(dims, NULL);
        int64_t *haloTags;
        haloTags = (int64_t *)malloc(sizeof(int64_t)*dims[0]);
        haloTags_set.read(haloTags, PredType::NATIVE_INT64, haloTags_space);

        rank = nodeMasses_space.getSimpleExtentDims(dims, NULL);
        double *nodeMasses;
        nodeMasses = (double *)malloc(sizeof(double)*dims[0]);
        nodeMasses_set.read(nodeMasses, PredType::NATIVE_DOUBLE, nodeMasses_space);

        rank = timesteps_space.getSimpleExtentDims(dims, NULL);
        int64_t *timesteps;
        timesteps = (int64_t *)malloc(sizeof(int64_t)*dims[0]);
        timesteps_set.read(timesteps, PredType::NATIVE_INT64, timesteps_space);

        rank = firstNodes_space.getSimpleExtentDims(dims, NULL);
        int64_t *firstNodes;
        firstNodes = (int64_t *)malloc(sizeof(int64_t)*dims[0]);
        firstNodes_set.read(firstNodes, PredType::NATIVE_INT64, firstNodes_space);

        rank = numberOfNodes_space.getSimpleExtentDims(dims, NULL);
        int64_t *numberOfNodes;
        numberOfNodes = (int64_t *)malloc(sizeof(int64_t)*dims[0]);
        numberOfNodes_set.read(numberOfNodes, PredType::NATIVE_INT64, numberOfNodes_space);

        // Import sodproperties, make dictionaries

        for(int i=0; i<dims[0]; i++)
        {
            // extract next tree from current hdf5
            int64_t first_ind = firstNodes[i];
            int64_t last_ind = firstNodes[i] + numberOfNodes[i];
    
            vector<int64_t> nodeIndex_this, descendentIndex_this, haloTags_this, timesteps_this;
            vector<double> nodeMasses_this;

            for(int jj=first_ind; jj<last_ind; jj++)
            {
                // get all tree members
                nodeIndex_this.push_back(nodeIndex[jj]);
                descendentIndex_this.push_back(descendentIndex[jj]);
                haloTags_this.push_back(haloTags[jj]);
                nodeMasses_this.push_back(nodeMasses[jj]);
                timesteps_this.push_back(timesteps[jj]);
                if(nodeMasses[jj] < min_z0_count_threshold)
                {
                    cout << nodeIndex[jj] << " " << descendentIndex[jj] << " " << haloTags[jj] << " " << nodeMasses[jj] << " " << timesteps[jj] << endl;
                }
            }

            set<int64_t> steps_present_set(timesteps_this.begin(), timesteps_this.end());
            vector<int> steps_present;
            for(set<int64_t>::iterator iter=steps_present_set.begin(); iter!=steps_present_set.end(); ++iter)
            {
                steps_present.push_back(*iter);
            }

            // -------- build backbone --------
            vector<int64_t> backbone_tags;
            vector<int> backbone_timesteps;
            int64_t backbone_member = -1;

            // find backbone root  at z=0
            for(int j=0; j<timesteps_this.size(); j++)
                if(timesteps_this[j] == final_timestep) 
                {
                     if(nodeMasses_this[j] >= min_z0_count_threshold*particlemass)
                     {
                         backbone_member = nodeIndex_this[j];
                         backbone_tags.push_back(haloTags_this[j]);

                         backbone_timesteps.push_back(timesteps_this[j]);
                         break;
                     }
                     else backbone_member = -2;
                }

            if(backbone_member == -1)
            {
                cout << "root not found" << endl;
                return 1;
            }

            // populate backbone toward higher redshift
            vector<int64_t> progenitor_indices_this_step;
            vector<int64_t> progenitor_tags_this_step;
            vector<double> progenitor_masses_this_step; 
            if(backbone_member > 0)
            {
                for(int j=steps_present.size()-2; j>=0; j--)
                {
                    int step = steps_present[j];
                    for(int k=0; k<timesteps_this.size(); k++)
                    {
                        if(timesteps_this[k] == step && descendentIndex_this[k] == backbone_member)
                        {
                            progenitor_indices_this_step.push_back(nodeIndex_this[k]);
                            progenitor_tags_this_step.push_back(haloTags_this[k]);
                            progenitor_masses_this_step.push_back(nodeMasses_this[k]);
                        }  
                    }
                    if(progenitor_masses_this_step.size() > 0)
                    {
                        double maxmass = 0.;
                        int maxmass_ind = -1;
                        for(int k=0; k<progenitor_masses_this_step.size(); k++)
                        {
                            if(progenitor_masses_this_step[k] > maxmass)
                            {
                                maxmass = progenitor_masses_this_step[k];
                                maxmass_ind = k;
                            }
                        }
                        backbone_member = progenitor_indices_this_step[maxmass_ind];
                        backbone_tags.push_back(progenitor_tags_this_step[maxmass_ind]);
                        backbone_timesteps.push_back(step);
                    }

                    else if(progenitor_masses_this_step.size() == 0) break;
                    progenitor_masses_this_step.clear();
                    progenitor_indices_this_step.clear();
                    progenitor_tags_this_step.clear();
                }

                for(int j=0; j<backbone_timesteps.size(); j++)
                {
                    root_index_vec.push_back(backbone_tags[0]);
                    step_vec.push_back(backbone_timesteps[j]);
                    hid_vec.push_back(backbone_tags[j]);
                
                    keykey_backbone_halos = make_pair(backbone_timesteps[j], abs(backbone_tags[j]));
                    backbone_halos_map.insert(make_pair(keykey_backbone_halos, true));
                }
            }
        }
        free(nodeIndex);
        free(descendentIndex); 
        free(haloTags);
        free(nodeMasses);
        free(timesteps);
        free(firstNodes);
        free(numberOfNodes);
    }

    map< pair<int, int64_t>, pair<double, double> > step_map;
    map< pair<int, int64_t>, double > step_map_r200;
    pair<int, int64_t> keykey;
    pair<double, double> valval;

    // SOD PROPERTIES FILE
    for(int step_ind=0; step_ind<count_steps; step_ind++)
    {
        int this_step = catalog_steps[step_ind];

        cout << endl << "loading sod properties step " << this_step << endl;

        tmp_string << "/projects/DarkUniverse_esp/heitmann/OuterRim/M000/L360/HACC001/analysis/Halos/M200/STEP" << this_step << "/m000-" << this_step << ".sodproperties";
        string haloproperties_gio_file = tmp_string.str();
        tmp_string.str("");

        int NR_H_properties;

        gio::GenericIO haloproperties_reader(MPI_COMM_SELF, haloproperties_gio_file);
        haloproperties_reader.openAndReadHeader(gio::GenericIO::MismatchRedistribute);
        int nRanks = haloproperties_reader.readNRanks();

        size_t MaxNElem_H_properties = 0;
        for (int j=0; j<nRanks; ++j) {
            size_t current_size = haloproperties_reader.readNumElems(j);
            MaxNElem_H_properties = current_size > MaxNElem_H_properties ? current_size : MaxNElem_H_properties;
        }

        int64_t *fof_halo_tag = new int64_t[ MaxNElem_H_properties ];
        int64_t *fof_halo_count = new int64_t[ MaxNElem_H_properties ];
        int64_t *sod_halo_count = new int64_t[ MaxNElem_H_properties ];
        float *sod_halo_radius = new float[ MaxNElem_H_properties ];
        float *sod_halo_mass = new float[ MaxNElem_H_properties ];
        float *sod_halo_cdelta = new float[ MaxNElem_H_properties ];
        float *sod_halo_cdelta_error = new float[ MaxNElem_H_properties ];
        float *sod_halo_c_acc_mass = new float[ MaxNElem_H_properties ];

        cout << "here" << endl;

        haloproperties_reader.addVariable( "fof_halo_tag", fof_halo_tag, true );
        haloproperties_reader.addVariable( "fof_halo_count", fof_halo_count, true );
        haloproperties_reader.addVariable( "sod_halo_count", sod_halo_count, true );
        haloproperties_reader.addVariable( "sod_halo_radius", sod_halo_radius, true );
        haloproperties_reader.addVariable( "sod_halo_mass", sod_halo_mass, true );
        haloproperties_reader.addVariable( "sod_halo_cdelta", sod_halo_cdelta, true );
        haloproperties_reader.addVariable( "sod_halo_cdelta_error", sod_halo_cdelta_error, true );
        haloproperties_reader.addVariable( "sod_halo_c_acc_mass", sod_halo_c_acc_mass, true );

        bool foundflag = false;

        int halos_found = 0;
        cout << "here" << endl;
        for(int j=0; j<nRanks; ++j) {
            size_t current_size = haloproperties_reader.readNumElems(j);
            cout << "Reading:" << current_size << endl;
            haloproperties_reader.readData(j);
            int nhalos = current_size; 
            for(int halo=0; halo<nhalos; halo++)
            {
                if(sod_halo_count[halo] >= min_sod_count)
                {
                    keykey = make_pair(this_step, fof_halo_tag[halo]);
                    valval = make_pair(sod_halo_count[halo], sod_halo_cdelta[halo]);
                    step_map.insert(make_pair(keykey, valval));
                    step_map_r200.insert(make_pair(keykey, sod_halo_radius[halo]));
                    halos_found ++;
                }
            }
        }

        haloproperties_reader.close();


        // SOD PROPERTY BINS FILES FOR CIRCULAR VELOCITY DEFINITIONS
        pair<int, int64_t> keykey_sod_property_bins;
        map< pair<int, int64_t>, pair<double, double> > backbone_halo_vcirc_map;
        map< pair<int, int64_t>, double > backbone_halo_rsmass_map;

        tmp_string << "/projects/DarkUniverse_esp/heitmann/OuterRim/M000/L360/HACC001/analysis/Halos/M200/STEP" << this_step << "/m000-" << this_step << ".sodpropertybins";
        string halopropertybins_gio_file = tmp_string.str();
        tmp_string.str("");

        int NR_H_propertybins;

        gio::GenericIO propertybins_reader(MPI_COMM_SELF, halopropertybins_gio_file);
        propertybins_reader.openAndReadHeader(gio::GenericIO::MismatchRedistribute);
        int nRanks_bins = propertybins_reader.readNRanks();

        size_t MaxNElem_H_bins = 0;
        for (int j=0; j<nRanks_bins; ++j) {
            size_t current_size = propertybins_reader.readNumElems(j);
            MaxNElem_H_bins = current_size > MaxNElem_H_bins ? current_size : MaxNElem_H_bins;
        }

        int64_t *fof_halo_bin_tag = new int64_t[ MaxNElem_H_bins ];
        int32_t *sod_halo_bin = new int32_t[ MaxNElem_H_bins ];
        int32_t *sod_halo_bin_count = new int32_t[ MaxNElem_H_bins ];
        float *sod_halo_bin_mass = new float[ MaxNElem_H_bins ];
        float *sod_halo_bin_radius = new float[ MaxNElem_H_bins ];

        propertybins_reader.addVariable( "fof_halo_bin_tag", fof_halo_bin_tag, true );
        propertybins_reader.addVariable( "sod_halo_bin", sod_halo_bin, true );
        propertybins_reader.addVariable( "sod_halo_bin_count", sod_halo_bin_count, true );
        propertybins_reader.addVariable( "sod_halo_bin_mass", sod_halo_bin_mass, true );
        propertybins_reader.addVariable( "sod_halo_bin_radius", sod_halo_bin_radius, true );

        for(int j=0; j<nRanks_bins; ++j) {
            size_t current_size = propertybins_reader.readNumElems(j);
            cout << "Reading:" << current_size << endl;
            propertybins_reader.readData(j);
            int nhalos = current_size / numbins;
            for(int halo=0; halo<nhalos; halo++)
            {
                int64_t getdata_IDs[numbins];
                int32_t getdata_halo_bin[numbins];
                float getdata_halo_bin_count[numbins];
                float getdata_halo_bin_radius[numbins];

                for(int bin=0; bin<numbins; bin++)
                {
                    getdata_IDs[bin] = fof_halo_bin_tag[halo*numbins + bin];
                    getdata_halo_bin[bin] = sod_halo_bin[halo*numbins + bin];
                    getdata_halo_bin_count[bin] = sod_halo_bin_count[halo*numbins + bin];
                    getdata_halo_bin_radius[bin] = sod_halo_bin_radius[halo*numbins + bin];
                }
      
                keykey_sod_property_bins=make_pair(this_step, abs(getdata_IDs[0]));
                if(backbone_halos_map.count(keykey_sod_property_bins) > 0)
                    if(backbone_halos_map[keykey_sod_property_bins] == true)
                    {
                        double this_halo_vpeak, this_halo_vr200, this_halo_rs_encmass;

                        
                        float sod_m200_count, c200;
                        float sod_r200; // fit only shells fully inside r200
                        if(step_map.count(keykey_sod_property_bins) > 0)
                        {
                            sod_m200_count = step_map[keykey_sod_property_bins].first;
                            sod_r200 = step_map_r200[keykey_sod_property_bins];
                            c200 = step_map[keykey_sod_property_bins].second;
                        }
                        else
                        {
                            cout <<  "ERROR halo not found" << endl;
                            exit(1);
                        }

                        float lowest_r = exp(log(getdata_halo_bin_radius[0]) - (log(getdata_halo_bin_radius[1])-log(getdata_halo_bin_radius[0]))); // inner radius of first shell - not given, calculated from width of other shells

                        double enclosedmass_masses[numbins]; // enclosed mass
                        double enclosedmass_radii[numbins];  // radius
                        double vcirc_profile[numbins];

                        for(int i=0; i<numbins; i++) // fill radii
                            enclosedmass_radii[i] = static_cast<double>(getdata_halo_bin_radius[i]);

                        for(int i=0; i<numbins; i++) // fill masses
                        {
                            float encmass_sum = 0;
                            for(int j=0; j<=i; j++)
                                encmass_sum += getdata_halo_bin_count[j];
                            enclosedmass_masses[i] = static_cast<double>(encmass_sum);
                        }

                        for(int i=1; i<numbins; i++)
                        {
                            if(enclosedmass_masses[i] <= enclosedmass_masses[i-1]) enclosedmass_masses[i] = enclosedmass_masses[i-1] + 1.e-10;
                        }

                        // linear interpolation of radius as a function of enclosed mass counts
                        gsl_interp *interp_enc_mass_r = gsl_interp_alloc(gsl_interp_linear, numbins);
                        gsl_interp_init(interp_enc_mass_r, enclosedmass_radii, enclosedmass_masses, numbins);
                        gsl_interp_accel *accelerator_enc_mass_r = gsl_interp_accel_alloc();

                        float m_of_r200 = gsl_interp_eval(interp_enc_mass_r, enclosedmass_radii, enclosedmass_masses, sod_r200, accelerator_enc_mass_r);
                        float diff = sod_m200_count - m_of_r200; // underestimation error - difference from linear interpolation, missing bin 0
                        for(int i=0; i<numbins; i++) // add back in the underestimation error
                            enclosedmass_masses[i] += diff;

                        // mass enclosed by r_s
                        if(c200 != -1 && sod_r200/c200 >= enclosedmass_radii[0] && sod_r200/c200 <= enclosedmass_radii[numbins-1])
                        {
                            gsl_interp *interp_encmass = gsl_interp_alloc(gsl_interp_linear, numbins);
                            gsl_interp_init(interp_encmass, enclosedmass_radii, enclosedmass_masses, numbins);
                            gsl_interp_accel *accelerator_encmass = gsl_interp_accel_alloc();
                            this_halo_rs_encmass = gsl_interp_eval(interp_encmass, enclosedmass_radii, enclosedmass_masses, sod_r200/c200, accelerator_encmass);
	
			    gsl_interp_free(interp_encmass);
			    gsl_interp_accel_free(accelerator_encmass);
                        }
                        else this_halo_rs_encmass = -1;

                        // Calculate circular velocity profile from enclosed masses
                        for(int i=0; i<numbins; i++)
                            vcirc_profile[i] = sqrt(G_gravity*particlemass*enclosedmass_masses[i]/enclosedmass_radii[i]/1000.); // 1000 to get radii in units of kpc instead of Mpc

                        // linear interpolation of circular velocity as a function of radius
                        this_halo_vr200 = -1;
                        gsl_interp *interp_vcirc = gsl_interp_alloc(gsl_interp_linear, numbins);
                        gsl_interp_init(interp_vcirc, enclosedmass_radii, vcirc_profile, numbins);
                        gsl_interp_accel *accelerator_vcirc = gsl_interp_accel_alloc();

                        this_halo_vr200 = gsl_interp_eval(interp_vcirc, enclosedmass_radii, vcirc_profile, sod_r200, accelerator_vcirc);
                        if(this_halo_vr200 < 0)
                        {
                            cout << "no vr200" << endl;
                            return 1;
                        }

                         // Find peak v_circ--no need for smoothing since this is a nice smooth function
                        this_halo_vpeak = 0;
                        int this_halo_vpeak_location = -1;
                        for(int i=0; i<numbins; i++)
                            if(vcirc_profile[i] > this_halo_vpeak)
                            {
                                this_halo_vpeak = vcirc_profile[i];
                                this_halo_vpeak_location = i;
                            }
                        if(this_halo_vpeak_location < 0)
                        {
                            cout << "no max found" << endl;
                            exit(1);
                        }

			gsl_interp_free(interp_enc_mass_r);
			gsl_interp_free(interp_vcirc);

			gsl_interp_accel_free(accelerator_enc_mass_r);
			gsl_interp_accel_free(accelerator_vcirc);

                        backbone_halo_vcirc_map.insert(make_pair(keykey_sod_property_bins, make_pair(this_halo_vpeak, this_halo_vr200)));
                        backbone_halo_rsmass_map.insert(make_pair(keykey_sod_property_bins, this_halo_rs_encmass));
                    }
              }
        }
	propertybins_reader.close();

        pair<int, int64_t> key_to_search;
        double found_rs_mass;
        pair<double, double> found_m200_c200;
        pair<double, double> found_vpeak_vr200;
        bool found_rel;
        for(int i=0; i<hid_vec.size(); i++)
        {
            if(step_vec[i] == this_step)
            {
                key_to_search = make_pair(step_vec[i], abs(hid_vec[i]));
                if(step_map.count(key_to_search) > 0)
                {
                    found_m200_c200 = step_map[key_to_search];
                    
                    if(backbone_halo_vcirc_map.count(key_to_search) > 0)
                    {
                        found_vpeak_vr200 = backbone_halo_vcirc_map[key_to_search];
                        found_rs_mass = backbone_halo_rsmass_map[key_to_search];
                    }
                    else found_vpeak_vr200 = make_pair(-1.,-1.);
                    lightcone_c << root_index_vec[i] << "\t" << step_vec[i] << "\t" << abs(hid_vec[i]) << "\t" << found_m200_c200.first << "\t" << found_m200_c200.second << "\t" << found_vpeak_vr200.first << "\t" << found_vpeak_vr200.second << "\t" << found_rs_mass << endl;
                }
            }
        }

        // Clear Maps
        backbone_halo_vcirc_map.clear();
        backbone_halo_rsmass_map.clear();
        step_map.clear();
        step_map_r200.clear();

        delete[] fof_halo_tag;
        delete[] fof_halo_count;
        delete[] sod_halo_count;
        delete[] sod_halo_radius;
        delete[] sod_halo_mass;
        delete[] sod_halo_cdelta;
        delete[] sod_halo_cdelta_error;
        delete[] sod_halo_c_acc_mass;

        delete[] fof_halo_bin_tag;
        delete[] sod_halo_bin;
        delete[] sod_halo_bin_count;
        delete[] sod_halo_bin_mass;
        delete[] sod_halo_bin_radius;

    }

    lightcone_c.close(); 

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
