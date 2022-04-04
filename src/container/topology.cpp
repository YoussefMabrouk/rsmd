/************************************************
 *                                              *
 *                rs@md                         *
 *    (reactive steps @ molecular dynamics )    *
 *                                              *
 ************************************************/
/* 
 Copyright 2020 Myra Biedermann
 Licensed under the Apache License, Version 2.0 
*/

#include "container/topology.hpp"
#include <iostream>
#include <math.h>      
using namespace std;

//
// get reaction record for specific molecule
//
const std::size_t& Topology::getReactionRecordMolecule(const std::size_t& oldmolid) 
{
    auto it = std::find_if(reactedMoleculeRecords.begin(), reactedMoleculeRecords.end(), [&oldmolid](const auto& rec){ return rec.first == oldmolid; });
    if( it == reactedMoleculeRecords.end() ) rsmdCRITICAL("couldn't find record for reacted molecule in topology: " << oldmolid);
    return it->second;
}

//
// get specific molecules
//
const Molecule& Topology::getMolecule(std::size_t molid) const
{
    // attention: assumes that molid is unique and thus returns first molecule that matches molid
    auto it = std::find_if( begin(), end(), [&molid](auto& m){ return molid == m.getID(); } );
    if( it == end() )   rsmdCRITICAL("couldn't find molecule in topology");
    return std::cref(*it);
}

std::vector<std::reference_wrapper<Molecule>> Topology::getMolecules(std::string molname)
{   
    // attention: returns all molecules that match molname
    std::vector<std::reference_wrapper<Molecule>> molReferences {};
    for(auto it = begin(); it != end(); it++ )
    {
        if( it->getName() == molname )  molReferences.emplace_back( *it );
    }
    return molReferences;
}


int Topology::heaviside(int i )
{
   if (i>0) 
   {
     return 1;
   }
   else 
   {
     return 0;
   }
}
int Topology::right( int n)
{ 
    return (n+1)*heaviside(CellNumbers[0]-1-n);
}
int Topology::left( int n)
{
    return (n-1)*heaviside(n) + (CellNumbers[0]-1)*heaviside(1-n);
}
int Topology::up( int n)
{
    return (n+1)*heaviside(CellNumbers[2]-1-n);
}
int Topology::down( int n)
{
    return (n-1)*heaviside(n) + (CellNumbers[2]-1)*heaviside(1-n);
}  
 
//3-d cell list
std::tuple<std::vector<std::vector<std::reference_wrapper<Molecule>>>, std::vector<std::vector<int>>> Topology::getCellList()
{
    std::vector<std::reference_wrapper<Molecule>> List {};
    std::vector<int> IndexList {};
    std::vector<std::vector<std::reference_wrapper<Molecule>>> CellList {};
    std::vector<std::vector<int>> CellNeighbourIndices {};
    Molecule molecule;
    int i, j, k, Index, NeighbourIndex;
    int n_x, n_y, n_z;
    for(i = 0 ; i < CellNumbers[0]*CellNumbers[1]*CellNumbers[2];i++)
    {
        CellList.emplace_back(List);
        CellNeighbourIndices.emplace_back(IndexList);
    }
    for(auto it = begin(); it != end(); it++ )
    {
        molecule = *it;
        n_x = floor((molecule[0].position(0)/dimensions[0]-floor(molecule[0].position(0)/dimensions[0]))*CellNumbers[0]);
        n_y = floor((molecule[0].position(1)/dimensions[1]-floor(molecule[0].position(1)/dimensions[1]))*CellNumbers[1]);
        n_z = floor((molecule[0].position(2)/dimensions[2]-floor(molecule[0].position(2)/dimensions[2]))*CellNumbers[2]); 
        Index =  n_x + n_y*CellNumbers[0] + n_z*CellNumbers[0]*CellNumbers[1];
        CellList[Index].emplace_back( *it );
    }
    for (k = 0; k<CellNumbers[2]; k++)
    {
        for (j = 0; j<CellNumbers[1]; j++)
        {
            for (i = 0; i<CellNumbers[0]; i++)
            {
                for ( auto n_x : {i, right(i), left(i)} )
                {
                    for ( auto n_y : {j, right(j), left(j)} ) 
                    {
                        for ( auto n_z : {k, up(k), down(k)} ) 
                        {
                            Index = i + j*CellNumbers[0] + k*CellNumbers[0]*CellNumbers[1];
                            NeighbourIndex = n_x +  n_y*CellNumbers[0] + n_z*CellNumbers[0]*CellNumbers[1];
                            CellNeighbourIndices[Index].emplace_back( NeighbourIndex );
                        }
                    }
                }
            }
        }
    }
    return {CellList, CellNeighbourIndices};
}


//
// get specific molecule and add it if not existing yet
//
Molecule& Topology::getAddMolecule(std::size_t molid, std::string molname)
{
    // attention: returns first molecule that matches molid (assumes that molid is unique)
    auto it = std::find_if( begin(), end(), [&molid](auto& m){ return molid == m.getID(); } );
    if( it == end() )
    {
        auto newIt = addMolecule( molid, molname );
        return std::ref(*newIt);
    }  
    else
        return std::ref(*it);
}



//
// remove specific molecule
//
void Topology::removeMolecule(Molecule& mol)
{
    data.erase( std::remove_if( begin(), end(), [&](auto& m){ return ( mol.getID() == m.getID() && mol.getName() == m.getName() ); } ), end() );
}

void Topology::removeMolecule(std::size_t molid)
{
    data.erase( std::remove_if( begin(), end(), [&](auto& m){ return molid == m.getID(); } ), end() );
}

//
// check if specific molecule exists in topology
//
bool Topology::containsMolecule(const Molecule& mol) const
{
    const auto tmp = std::find_if( begin(), end(), [&](auto& m){ return ( mol.getID() == m.getID() && mol.getName() == m.getName() ); } );
    if( tmp != end() )
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Topology::containsMolecule(const std::size_t& molid) const
{
    const auto tmp = std::find_if( begin(), end(), [&](auto& m){ return molid == m.getID(); } );
    if( tmp != end() )
    {
        return true;
    }
    else
    {
        return false;
    }

}

//
// get moleculetypes
//
std::vector<std::string> Topology::getMoleculetypes() const
{
    std::vector<std::string> moleculetypes;
    for( const auto& m: data )
    {
        auto it = std::find_if( moleculetypes.begin(), moleculetypes.end(), [&m](auto mt){ return mt == m.getName(); } );
        if( it == moleculetypes.end() )    moleculetypes.push_back( m.getName() );
    }
    return moleculetypes;

}

//
// sort topology, i.e. rearrange and renumber everything (molecules + atoms)
//
void Topology::sort()
{
    // clear atomic reaction records
    reactedAtomRecords.clear();

    // sort (according to name) and renumber molecules
    // then renumber atoms accordingly
    // note:  use stable_sort instead of sort to retain order of equal elements!
    std::stable_sort( begin(), end(), [](auto& lhs, auto& rhs){ return lhs.getName() < rhs.getName(); });
    
    std::size_t counterMolecules = 0;
    std::size_t counterAtoms = 0;
    
    for( auto& m: data )
    {   
        //cout<<"aaah" << m << endl;
        // renumber molecules
        ++ counterMolecules;
        // check if this is a newly reacted molecule
        bool isReactedMolecule = false;
        auto search = std::find_if(std::begin(reactedMoleculeRecords), std::end(reactedMoleculeRecords), [&m](const auto& record){return m.getID() == record.first; });
        if( search != std::end(reactedMoleculeRecords) )
        {
            isReactedMolecule = true;
            search->second = counterMolecules;
        }
        // reset ID
        #ifndef NDEBUG
        if( m.getID() != counterMolecules ){ rsmdDEBUG("note: resetting ID of " << m << " to " << counterMolecules); }
        #endif
        m.setID(counterMolecules);
        // renumber atoms in molecule
        for( auto& a: m )
        {
            ++ counterAtoms;
            // record ID changes if reactedMolecule
            if( isReactedMolecule ) reactedAtomRecords.push_back(std::make_pair(a.id, counterAtoms));
            // update ID
            #ifndef NDEBUG
            if( a.id != counterAtoms ){ rsmdDEBUG("note: resetting ID of " << a << " to " << counterAtoms ); }
            #endif
            a.id = counterAtoms;
        }
    }
}

//
// repair (new) molecules that are broken across periodic boundaries
//
void Topology::repairMoleculePBC(Molecule& molecule)
{
    rsmdDEBUG( "repairing molecule, in case it is broken across periodic boundaries: " << molecule );
    Atom& referenceAtom = molecule.front();
    for(auto& atom: molecule)
    {
        rsmdDEBUG( "   before: " << atom );
        REALVEC distance = (atom.position - referenceAtom.position);
        for( std::size_t i=0; i<3; ++i )
        {
            atom.position[i] -= static_cast<int>( distance[i] / (0.5 * dimensions[i]) ) * dimensions[i];
        }
        rsmdDEBUG( "   after: " << atom );
    }
}
