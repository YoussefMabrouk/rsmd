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

#include "container/universe.hpp"
#include "container/topology.hpp"
#include <iostream>
#include <cmath>
#include <math.h>
using namespace std;
//
// initial setup of the universe
//
void Universe::setup(const Parameters& parameters)     
{
    switch( parameters.getEngineType() )
    {
        case ENGINE::GROMACS:   
            topologyParser = std::make_unique<TopologyParserGMX>();
            assert(topologyParser);

            unitSystem = std::make_unique<UnitSystem>("nm", "ps", "kJ/mol", "K");
            assert(unitSystem);

            break;

        case ENGINE::NONE:
            rsmdCRITICAL( "md engine is set to none" );
            break;
    }

    // read reaction templates from files
    auto reactionFiles = parameters.getOption("reaction.file").as<std::vector<std::string>>();
    rsmdLOG("... reading reaction templates ... ");
    ReactionParser reactionParser {};
    for( const auto& file: reactionFiles )
    {
        auto reaction = reactionParser.read( file );
        
        // some verbose printing: 
        rsmdLOG( "... from file '" << file << "': ")
        rsmdLOG( reaction );

        rsmdLOG( "... checking for consistency in provided input for reaction '" << reaction.getName() << "' ...");
        // check that reaction template contains required input 
        // for the chosen simulation algorithm
        switch( parameters.getSimulationAlgorithm() )
        {
            case SIMALGORITHM::MC:
                // --> check for energy != 0
                if( reaction.getReactionEnergy() == 0 )
                    rsmdWARNING( "    reaction energy == 0, are you sure that is correct?" );
                break;

            case SIMALGORITHM::RATE:
                // --> check for reaction rate
                if( reaction.getRate().size() == 0 )
                    rsmdWARNING( "    no reaction rate input, are you sure that is correct?" );
                break;
        }
        // check for consistency within reactants/products/criterions
        reaction.consistencyCheck();
        rsmdLOG( "... consistency check done. everything seems fine.");
        
        reactionTemplates.emplace_back(reaction);
    }
}


// 
// update topologies
//
void Universe::update(const std::size_t& cycle) 
{
    topologyOld.clear();
    topologyNew.clear();
    topologyRelaxed.clear();
    topologyParser->read(topologyOld, cycle);
    topologyOld.clearReactionRecords();
    topologyNew = topologyOld;
}


//
// write (new) topology to file 
//
void Universe::write(const std::size_t& cycle)
{
    topologyNew.sort();
    topologyParser->write(topologyNew, cycle);
}


//
// read relaxed configuration from file
//
void Universe::readRelaxed(const std::size_t& cycle)
{
    topologyRelaxed.clear();
    topologyParser->readRelaxed(topologyRelaxed, cycle);
}


//
// check given reaction candidate (in topologyRelaxed) for 
// 'physical meaningfulness' after relaxation, 
// i.e. how much the corresponding atoms moved
//
void Universe::checkMovement(const ReactionCandidate& candidate)
{
    // first: compute typical length in system against which to check
    REAL volume = topologyNew.getDimensions()[0] * topologyNew.getDimensions()[1] * topologyNew.getDimensions()[2];
    // REAL typicalDistance = std::sqrt( (3.0 * volume) / (4.0 * M_PI * topologyNew.getNAtoms()) );
    REAL typicalDistance = std::cbrt( (3.0 * volume) / (4.0 * M_PI * topologyNew.getNAtoms()) );

    for( auto& molecule: candidate.getProducts() )
    {
        // get same molecule in topologyRelaxed
        std::size_t newMolID = topologyNew.getReactionRecordMolecule(molecule.getID());
        auto& newMolecule = topologyRelaxed.getMolecule(newMolID);

        // go through molecule and compute movement of atoms
        auto atomBefore = molecule.begin();
        auto atomAfter  = newMolecule.begin();
        while( atomBefore != molecule.end() || atomAfter != newMolecule.end() )
        {
            auto distance = enhance::distance(*atomBefore, *atomAfter, topologyNew.getDimensions());

            if( distance > 3 * typicalDistance )
            {
                rsmdWARNING( std::setprecision(3) << "... atom " << atomAfter->name << " " << atomAfter->id << " of molecule " << newMolecule.getName() << " " << newMolecule.getID() 
                        << " moved more than three times the typical distance: " << distance << ' ' << unitSystem->length << " ( > 3 * " << typicalDistance << ' ' << unitSystem->length << ")");
            }
            else if( distance > 2 * typicalDistance )
            {
                rsmdWARNING( std::setprecision(3) << "... atom " << atomAfter->name << " " << atomAfter->id << " of molecule " << newMolecule.getName() << " " << newMolecule.getID() 
                        << " moved more than twice the typical distance: " << distance << ' ' << unitSystem->length << " ( > 2 * " << typicalDistance << ' ' << unitSystem->length << ")");
            }
            else
            {
                rsmdDEBUG( "... atom " << atomAfter->name << " " << atomAfter->id << " of molecule " << newMolecule.getName() << " " << newMolecule.getID() 
                        << " moved: " << distance << ' ' << unitSystem->length);
            }
            ++ atomBefore;
            ++ atomAfter;
        }
    }
}


//
// check if a candidate is still available
//
bool Universe::isAvailable( const ReactionCandidate& candidate )
{
    bool reactantsAreAvailable = true;
    for( const auto& reactant: candidate.getReactants() )
    {
        if( ! topologyNew.containsMolecule(reactant) )
        {
            rsmdDEBUG( "couldn't find molecule " << reactant.getName() << " " << reactant.getID() << " in topology" );
            reactantsAreAvailable = false;
            break;
        }
    }
    return reactantsAreAvailable;
}

void Universe::makeMoleculeWhole(Molecule& molecule, const REALVEC& dimensions)
{
    rsmdLOG( "... repairing molecule in case it is broken across periodic boundaries: " << molecule );
    Atom& referenceAtom = molecule.front();
    for(auto& atom: molecule)
    {   
        REALVEC before = atom.position;
        REALVEC distance = (atom.position - referenceAtom.position);
        bool moved = false;
        for( std::size_t i=0; i<3; ++i )
        {
            atom.position[i] -= static_cast<int>( distance[i] / (0.5 * dimensions[i]) ) * dimensions[i];
            if( static_cast<int>(distance[i] / (0.5*dimensions[i]) ) != 0 )    moved = true;
        }
        REALVEC after = atom.position;
        if( moved )
        {
            rsmdLOG( "    before: " << before );
            rsmdLOG( "    after: " << after );
        }
    }
}

//
// react a given candidate
// (checks for whether the molecules are still available need to happen before!)
//
void Universe::react(ReactionCandidate& candidate)
{
    rsmdDEBUG( "performing reaction for candidate " << candidate.shortInfo() );
   
    // reactant --> product translation 
    candidate.applyTransitions();
    // make products whole
    for(auto& product: candidate.getProducts())
    {
        makeMoleculeWhole(product, topologyNew.getDimensions());
    }
    // apply translational movements of product atoms
    candidate.applyTranslations();

    // apply changes to topology
    auto highestMolID = std::max_element( std::begin(topologyNew), std::end(topologyNew), [](const auto& mol1, const auto& mol2){ return mol1.getID() < mol2.getID(); } )->getID(); 
    for( const auto& reactant: candidate.getReactants() )
    {
        topologyNew.removeMolecule( reactant.getID() );    
    }
    for( auto& product: candidate.getProducts() )
    {
        product.setID( ++highestMolID );
        auto molecule __attribute__((unused)) = topologyNew.addMolecule( product );
        topologyNew.addReactionRecord( highestMolID );
        // topologyNew.repairMoleculePBC( *molecule );
        rsmdDEBUG( "new molecule " << molecule->getName() << " got ID " << molecule->getID() );
    }
}

//cell list 
std::tuple<std::vector<std::reference_wrapper<Molecule>>, std::vector<int>> Universe::CellNeighbours(int CellIndex , std::string molname)
{   
    std::vector<std::reference_wrapper<Molecule>> molReferences {};
    std::vector<int> molCells {};
    Molecule molecule;    
    int i, j, Index;

    for (i= 0 ; i < CellNeighbourIndices[CellIndex].size(); i++)
    {
      Index = CellNeighbourIndices[CellIndex][i];
      for(j = 0 ; j < CellList[Index].size(); j++)
      {
        molecule = CellList[Index][j];
        if( molecule.getName() == molname )  
        {
            molReferences.emplace_back( CellList[Index][j] );
            molCells.emplace_back( Index );
        }
      } 
    }
    return {molReferences, molCells};
}

std::vector<std::reference_wrapper<Molecule>> Universe::Cell(int CellIndex , std::string molname)
{   
    std::vector<std::reference_wrapper<Molecule>> molReferences {};
    Molecule molecule;   
    int j;
    
    for(j = 0 ; j < CellList[CellIndex].size(); j++)
    {
        molecule = CellList[CellIndex][j];
        if( molecule.getName() == molname )  molReferences.emplace_back( CellList[CellIndex][j] );
    }
    return molReferences;
}

std::vector<ReactionCandidate> Universe::CellSearchReactionCandidates()
{
    int i, CellIndex;
    std::vector<ReactionCandidate> reactionCandidates {};
    std::vector<double> reactionRates {};
    auto [x, y] = topologyOld.getCellList();
    CellList = x;
    CellNeighbourIndices = y;
    for(CellIndex = 0; CellIndex < CellList.size();CellIndex++)
    {
        for( auto& candidate: CellReactionCandidates ( CellIndex ))
        {
            reactionCandidates.push_back (candidate);
        }
    }
    //for(i = 0 ; i < reactionCandidates.size();i++)
    //{
    //    reactionRates.push_back(reactionCandidates[i].getCurrentReactionRateValue());
    //}
    enhance::weighted_shuffle(reactionCandidates.begin(), reactionCandidates.end(), reactionRates.begin(), reactionRates.end());
    
    return reactionCandidates;
}

std::vector<ReactionCandidate> Universe::CellReactionCandidates(int CellIndex)
{
    // search for possible reaction candidates and return them if they match all criteria
    std::vector<ReactionCandidate> reactionCandidates {};
    std::vector<std::reference_wrapper<Molecule>> reactants1;
    std::vector<std::reference_wrapper<Molecule>> reactants2;
    std::vector<std::reference_wrapper<Molecule>> reactants3;
    std::vector<std::reference_wrapper<Molecule>> reactants4;
    std::vector<int> CellIndex2, CellIndex3, CellIndex4;
    Molecule reactant1;
    Molecule reactant2;
    Molecule reactant3;
    Molecule reactant4;
    int i, j, k, l, cellindex1, cellindex2, cellindex3, cellindex4;
    
    for( auto& reactionTemplate: reactionTemplates )
    {
        if( reactionTemplate.getReactants().size() == 2 )
        {            
            reactants1 = Cell(CellIndex, reactionTemplate.getReactants()[0].getName() );
            for(i = 0 ; i < reactants1.size();i++)
            {
              reactant1 = reactants1[i];
              reactionCandidates.push_back( reactionTemplate );
              reactionCandidates.back().updateReactant( 0, reactant1 );
              rsmdDEBUG( "checking reaction candidate: " << reactant1.getName() << ", " << reactant1.getID() );
              if( reactionCandidates.back().valid(topologyOld.getDimensions(), 0) )
              {
                  reactionCandidates.pop_back();
                  auto [reactants2, CellIndex2] = CellNeighbours(CellIndex, reactionTemplate.getReactants()[1].getName() );
                  for(j = 0 ; j < reactants2.size();j++)
                  {
                      reactant2 = reactants2[j];
                      if( reactant1.getID() == reactant2.getID() ) continue;
                      if( reactant1.getName() == reactant2.getName() && reactant1.getID() > reactant2.getID() ) continue;
                      rsmdDEBUG( "checking reaction candidate: " << reactant2.getName() << ", " << reactant2.getID() );
                      reactionCandidates.push_back( reactionTemplate );
                      reactionCandidates.back().updateReactant( 0, reactant1 );
                      reactionCandidates.back().updateReactant( 1, reactant2 );
                      if( ! reactionCandidates.back().valid(topologyOld.getDimensions(), 1) ) reactionCandidates.pop_back();
                  }
              }
              else
              {
                reactionCandidates.pop_back();  
              }
            }
        }         
        else if( reactionTemplate.getReactants().size() == 3 )
        {
            reactants1 = Cell(CellIndex, reactionTemplate.getReactants()[0].getName() );
            for(i = 0 ; i < reactants1.size();i++)
            {
              reactant1 = reactants1[i];
              reactionCandidates.push_back( reactionTemplate );
              reactionCandidates.back().updateReactant( 0, reactant1 );
              rsmdDEBUG( "checking reaction candidate: " << reactant1.getName() << ", " << reactant1.getID() );
              if ( reactionCandidates.back().valid(topologyOld.getDimensions(), 0))
              {
                  reactionCandidates.pop_back();
                  auto [reactants2, CellIndex2] = CellNeighbours(CellIndex, reactionTemplate.getReactants()[1].getName() );
                  for(j = 0 ; j < reactants2.size();j++)
                  {
                      reactant2 = reactants2[j];
                      if( reactant1.getID() == reactant2.getID() ) continue;
                      if( reactant1.getName() == reactant2.getName() && reactant1.getID() > reactant2.getID() ) continue;
                      rsmdDEBUG( "checking reaction candidate: " << reactant2.getName() << ", " << reactant2.getID() );
                      reactionCandidates.push_back( reactionTemplate );
                      reactionCandidates.back().updateReactant( 0, reactant1 );
                      reactionCandidates.back().updateReactant( 1, reactant2 );
                      if( reactionCandidates.back().valid(topologyOld.getDimensions(), 1) )
                      {
                          reactionCandidates.pop_back();
                          auto [reactants3, CellIndex3] = CellNeighbours(CellIndex, reactionTemplate.getReactants()[2].getName() );
                          for(k = 0 ; k < reactants3.size();j++)
                          {
                             reactant3 = reactants3[k]
                             if( reactant1.getID() == reactant3.getID() | reactant2.getID() == reactant3.getID() ) continue;
                             if( reactant1.getName() == reactant3.getName() && reactant1.getID() > reactant3.getID() ) continue;
                             if( reactant2.getName() == reactant3.getName() && reactant2.getID() > reactant3.getID() ) continue;
                             reactionCandidates.push_back( reactionTemplate );
                             reactionCandidates.back().updateReactant( 0, reactant1 );
                             reactionCandidates.back().updateReactant( 1, reactant2 );                                             
                             reactionCandidates.back().updateReactant( 2, reactant3 );
                             if( ! reactionCandidates.back().valid(topologyOld.getDimensions(), 2) ) reactionCandidates.pop_back();
                          }
                       }
                       else
                       {
                         reactionCandidates.pop_back();
                       }                               
                  }
              }
              else
              {
                  reactionCandidates.pop_back();
              }
            }
         }        
        if( reactionTemplate.getReactants().size() == 4 )
        {
            reactants1 = Cell(CellIndex, reactionTemplate.getReactants()[0].getName() );
            for(i = 0 ; i < reactants1.size();i++)
            {
              reactant1 = reactants1[i];
              cellindex1 = CellIndex;
              reactionCandidates.push_back( reactionTemplate );
              reactionCandidates.back().updateReactant( 0, reactant1 );
              rsmdDEBUG( "checking reaction candidate: " << reactant1.getName() << ", " << reactant1.getID() );
              if ( reactionCandidates.back().valid(topologyOld.getDimensions(), 0))
              {
                  reactionCandidates.pop_back();
                  auto [reactants2, CellIndex2] = CellNeighbours(CellIndex, reactionTemplate.getReactants()[1].getName() );
                  for(j = 0 ; j < reactants2.size();j++)
                  {
                      reactant2 = reactants2[j];
                      cellindex2 = CellIndex2[j];
                      if( reactant1.getID() == reactant2.getID() ) continue;
                      if( reactant1.getName() == reactant2.getName() && reactant1.getID() > reactant2.getID() ) continue;   
                      if( reactant1.getName() == reactant2.getName() && cellindex1 > cellindex2 ) continue;
                      rsmdDEBUG( "checking reaction candidate: " << reactant2.getName() << ", " << reactant2.getID() );
                      reactionCandidates.push_back( reactionTemplate );
                      reactionCandidates.back().updateReactant( 0, reactant1 );
                      reactionCandidates.back().updateReactant( 1, reactant2 );
                      if( reactionCandidates.back().valid(topologyOld.getDimensions(), 1) )
                      {
                          reactionCandidates.pop_back();
                          auto [reactants3, CellIndex3] = CellNeighbours(CellIndex, reactionTemplate.getReactants()[2].getName() );
                          for (k = 0 ; k < reactants3.size();k++)
                          {
                              reactant3 = reactants3[k];
                              cellindex3 = CellIndex3[k];
                              if( reactant1.getID() == reactant3.getID() || reactant2.getID() == reactant3.getID() )  continue;
                              if( reactant1.getName() == reactant3.getName() && reactant1.getID() > reactant3.getID() ) continue;
                              if( reactant1.getName() == reactant3.getName() && cellindex1 > cellindex3 ) continue;
                              if( reactant2.getName() == reactant3.getName() && reactant2.getID() > reactant3.getID() ) continue;
                              if( reactant2.getName() == reactant3.getName() && cellindex2 > cellindex3 ) continue;
                              rsmdDEBUG( "checking reaction candidate: " << reactant3.getName() << ", " << reactant3.getID() );
                              reactionCandidates.push_back( reactionTemplate );
                              reactionCandidates.back().updateReactant( 0, reactant1 );
                              reactionCandidates.back().updateReactant( 1, reactant2 );                                             
                              reactionCandidates.back().updateReactant( 2, reactant3 );
                              if( reactionCandidates.back().valid(topologyOld.getDimensions(), 2) )
                              {
                                  reactionCandidates.pop_back(); 
                                  auto [reactants4, CellIndex4] = CellNeighbours(CellIndex, reactionTemplate.getReactants()[3].getName() );
                                  for (l = 0 ; l < reactants4.size();l++)
                                  {
                                      reactant4 = reactants4[l];
                                      cellindex4 = CellIndex4[l];
                                      if( reactant1.getID() == reactant4.getID() || reactant2.getID() == reactant4.getID() || reactant3.getID() == reactant4.getID() )  continue;
                                      if( reactant1.getName() == reactant4.getName() && reactant1.getID() > reactant4.getID() ) continue;
                                      if( reactant1.getName() == reactant4.getName() && cellindex1 > cellindex4 ) continue;
                                      if( reactant2.getName() == reactant4.getName() && reactant2.getID() > reactant4.getID() ) continue;
                                      if( reactant2.getName() == reactant4.getName() && cellindex2 > cellindex4 ) continue;
                                      if( reactant3.getName() == reactant4.getName() && reactant3.getID() > reactant4.getID() ) continue;
                                      if( reactant3.getName() == reactant4.getName() && cellindex3 > cellindex4 ) continue;
                                      rsmdDEBUG( "checking reaction candidate: " << reactant4.getName() << ", " << reactant4.getID() );
                                      reactionCandidates.push_back( reactionTemplate );
                                      reactionCandidates.back().updateReactant( 0, reactant1 );
                                      reactionCandidates.back().updateReactant( 1, reactant2 );
                                      reactionCandidates.back().updateReactant( 2, reactant3 );
                                      reactionCandidates.back().updateReactant( 3, reactant4 );
                                      if( ! reactionCandidates.back().valid(topologyOld.getDimensions(), 3) ) reactionCandidates.pop_back();
                                  }
                              }
                              else
                              {
                                  reactionCandidates.pop_back();
                              }
                          }    
                      }
                      else
                      {
                          reactionCandidates.pop_back();
                      }
                  }
              }
              else
              {
                  reactionCandidates.pop_back();
              }
            }
        }       
    }
    
    return reactionCandidates;
}

    
