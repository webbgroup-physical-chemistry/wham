#include "bw_hdf5.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

void Interact_H5::h5_init(bw_options option)
{
    options = option;

    std::stringstream ss;
    char outname[1024];
    sprintf(outname,"%s",options.hdf5file.c_str());
    ss << outname;
    ss >> filename;
    const H5std_string FILE_NAME(outname);
    try
    {
        Exception::dontPrint();
        file = new H5File(FILE_NAME, H5F_ACC_RDWR);
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

void Interact_H5::h5_write_dat(bw_datfile filedat)
{
    try
    {
        Exception::dontPrint();
        DataSpace *dataspace;
        DataSet *dataset;
        float * values;
        char name[1024];
        int ncol = filedat.dat[0].size();
        // For each column
        for (int i=0; i<ncol; i++)
        {
            int datsize = filedat.dat.size();
            values = new float [datsize];
            for (int j=0; j<datsize; j++)
            {
                values[j] = filedat.dat[j][i];
            }
            /* Write data */
            sprintf(name,"/Trajectories/Traj-%i/%s",filedat.experiment,options.datnames[i].c_str());
            /* 
               See if it exists already.  If it does, unlink it so we can replace it.
               This is set up as a while loop so that it can trivially be replaced with
               a file->move(src,dst) command to keep the old data.  I chose not to go
               that route, however, for fear of cluttering up the hdf5 file. Regardless, 
               the hdf5 file WILL continue to grow in size...
            */
            bool data_already_exists = true;
            int n=0;
            while (data_already_exists)
            {
                if (n>100){
                    std::cerr << "\nERROR: Found " << name << " " << " times.  Cowardly refusing to keep trying.\n" << std::endl;
                    std::exit(1);}
                try
                {
                    file->openDataSet(name);
                    file->unlink(name);
                }
                catch(...)
                {
                    data_already_exists = false;
                }
                n++;
            }
            hsize_t dim[1] = { datsize };
            dataspace = new DataSpace(1,dim);
            dataset = new DataSet(file->createDataSet(name,PredType::NATIVE_FLOAT,*dataspace));
            dataset->write(values,PredType::NATIVE_FLOAT);

            /* Attribute */
            StrType str_type(PredType::C_S1,1024);
            hsize_t attrsize = 1;
            char unit_attr[attrsize][1024];
            hsize_t attr_dim[1] = { attrsize };
            sprintf(unit_attr[0],"%s",options.datunits[i].c_str());
            DataSpace attr_dataspace(1,attr_dim);
            Attribute unit_attribute = dataset->createAttribute("Units",str_type,attr_dataspace);
            unit_attribute.write(str_type,unit_attr);
            std::cout << name << std::endl;
            
            dataset->close();
            delete dataset;
            delete dataspace;
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

int Interact_H5::h5_get_dataset(std::string path, h5_dat &data)
{
    try
    {
        Group grp = file->openGroup(path.c_str());
        // find out which frames this covered
        hsize_t att_dim[1] = { 2 };
        int attr_data[2];
        Attribute span = grp.openAttribute("First frame read, last frame read");
        span.read(PredType::NATIVE_INT,attr_data);
        data.span.push_back(attr_data[0]);
        data.span.push_back(attr_data[1]);
    }
    catch(...)
    {
        return -1;
    }
    char bin_name[1024],prob_name[1024];
    sprintf(bin_name,"%s/BinCounts",path.c_str());
    sprintf(prob_name,"%s/Probability",path.c_str());
    DataSet *dataset;
    DataSpace *dataspace;
    try
    {
        Exception::dontPrint();
        // read the bin counts
        dataset = new DataSet(file->openDataSet(bin_name));
        dataspace = new DataSpace(dataset->getSpace());
        int rank = dataspace->getSimpleExtentNdims();
        hsize_t dims[rank];
        int ndims = dataspace->getSimpleExtentDims(dims,NULL); 
        int bin_out[dims[0]];
        dataset->read(bin_out, PredType::NATIVE_INT);
        delete dataspace;
        delete dataset;
        // read the probabilities
        dataset = new DataSet(file->openDataSet(prob_name));
        dataspace = new DataSpace(dataset->getSpace());
        int prank = dataspace->getSimpleExtentNdims();
        hsize_t pdims[prank];
        int pndims = dataspace->getSimpleExtentDims(pdims,NULL);
        float prob_out[pdims[0]][pdims[1]];
        dataset->read(prob_out, PredType::NATIVE_FLOAT);
        delete dataspace;
        delete dataset; 
        // Make sure the sizes match
        if (pdims[0] != dims[0])
        {
            std::cerr << "\nERROR! " << bin_name << " and " << prob_name << " do not have the same number of bins!\n";
            std::exit(1);
        }
        /* 
           The probability of an individual bin is prob_out. But we 
           divide by the number of times that bin is visited to get 
           the probability of seeing a single frame in that bin.
        */
        for (int i=0;i<dims[0];i++)
        {
            if (bin_out[i] > 0)
            {
                // There are pdims[1] degrees of freedom, so we want the index of ndof - 1
                //data.bin_prob.push_back(prob_out[i][pdims[1]-1]/(float)bin_out[i]);
                /* Don't divide by the number of times the bin is visited yet */
                data.bin_prob.push_back(prob_out[i][pdims[1]-1]);
            }
            else
            {
                data.bin_prob.push_back(0);
            }
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return 0;
}

void Interact_H5::h5_bin_assignments(std::vector<bw_datfile> &list)
{
    /* Obtain the bin number for each frame */
    char name[1024];
    DataSet *dataset;
    DataSpace *dataspace;
    try
    {
        Exception::dontPrint();
        for (int i=0; i<(int)list.size(); i++)
        {
            sprintf(name,"/Trajectories/Traj-%i/Bin",list[i].experiment);
            //std::cout << name << std::endl;
            dataset = new DataSet(file->openDataSet(name));
            dataspace = new DataSpace(dataset->getSpace());
            int rank = dataspace->getSimpleExtentNdims();
            hsize_t dims[rank];
            int ndims = dataspace->getSimpleExtentDims(dims,NULL);
            int bin_out[dims[0]];
            dataset->read(bin_out,PredType::NATIVE_INT);
            delete dataspace;
            delete dataset;
            for (int n=0; n<dims[0]; n++)
            {
                list[i].bin.push_back(bin_out[n]);
            }
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

void Interact_H5::h5_write_prob(std::vector<h5_dat> probs)
{
    int nnn=0;
    try
    {
        Exception::dontPrint();
        DataSpace *dataspace;
        DataSet *dataset;
        char name[1024];
        // For each column
        int ncol = probs[0].avg.size();
        int nconv = probs.size();
        for (int i=0; i<ncol; i++)
        {
            float averages[nconv][2];
            for (int j=0; j<nconv; j++)
            {
                averages[j][0] = probs[j].avg[i];
                averages[j][1] = probs[j].stdev[i];
                //stdevs[j] = probs[j].stdev[i];
                //std::cout << i << " " << j << " " <<  options.datnames[i] << " " << averages[j][0] << "(" << averages[j][1] << ")" << "\n";
            }

            /* Write data */
            sprintf(name,"/Ensemble/%s",options.datnames[i].c_str());
            std::cout << "Writing " << name << std::endl;
            /*
             See if it exists already.  If it does, unlink it so we can replace it.
             This is set up as a while loop so that it can trivially be replaced with
             a file->move(src,dst) command to keep the old data.  I chose not to go
             that route, however, for fear of cluttering up the hdf5 file. Regardless,
             the hdf5 file WILL continue to grow in size...
             */
            bool data_already_exists = true;
            int n=0;
            while (data_already_exists)
            {
                if (n>100){
                    std::cerr << "\nERROR: Found " << name << " " << " times.  Cowardly refusing to keep trying.\n" << std::endl;
                    std::exit(1);}
                try
                {
                    file->openDataSet(name);
                    file->unlink(name);
                }
                catch(...)
                {
                    data_already_exists = false;
                }
                n++;
            }
            hsize_t dim[2] = { nconv, 2 };
            dataspace = new DataSpace(2,dim);
            dataset = new DataSet(file->createDataSet(name,PredType::NATIVE_FLOAT,*dataspace));
            dataset->write(averages,PredType::NATIVE_FLOAT);

            /* Attribute */
            StrType str_type(PredType::C_S1,1024);
            hsize_t attrsize = 2;
            char unit_attr[attrsize][1024];
            hsize_t attr_dim[1] = { attrsize };
            sprintf(unit_attr[0],"average(%s)",options.datunits[i].c_str());
            sprintf(unit_attr[1],"stdev(%s)",options.datunits[i].c_str());
            DataSpace attr_dataspace(1,attr_dim);
            Attribute unit_attribute = dataset->createAttribute("Units",str_type,attr_dataspace);
            unit_attribute.write(str_type,unit_attr);

            dataset->close();
            delete dataset;
            delete dataspace;
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;

}
