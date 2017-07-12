#include "InterfaceCellIdentifierImp.h"
#include <functional>
#include <algorithm>
#include <cmath>

#include <iostream>
using namespace std;

constexpr bool InterfaceCellIdentifierImp::IDENTIFIED;
constexpr bool InterfaceCellIdentifierImp::UNIDENTIFIED;


bool
InterfaceCellIdentifierImp::CellCondition::IsSame (
	const CellCondition& other
) const noexcept
{
	auto thisNNeighsIter = nNeigh2Other.begin();
	auto thisDiffIter = neighDiff.begin();
	auto otherNNeighsIter = other.nNeigh2Other.begin();
	auto otherDiffIter = other.neighDiff.begin();
	for ( ;	thisNNeighsIter != nNeigh2Other.end();
		++thisNNeighsIter, ++thisDiffIter,
		++otherNNeighsIter, ++otherDiffIter
	)
	{
		if ( *thisNNeighsIter != *otherNNeighsIter )
			return false;

		for ( const auto& thisToOuterDiff : *thisDiffIter )
			if ( !MatchAny(thisToOuterDiff, *otherDiffIter) )
				return false;
	}
	return true;
}


bool
InterfaceCellIdentifierImp::CellCondition::IsSimilar (
	const CellCondition& other
) const noexcept
{
	auto thisNNeighsIter = nNeigh2Other.begin();
	auto thisDiffIter = neighDiff.begin();
	auto otherNNeighsIter = other.nNeigh2Other.begin();
	auto otherDiffIter = other.neighDiff.begin();
	for ( ;	thisNNeighsIter != nNeigh2Other.end();
		++thisNNeighsIter, ++thisDiffIter,
		++otherNNeighsIter, ++otherDiffIter
	)
	{
		if ( otherDiffIter->empty() )
			continue;

		if ( *thisNNeighsIter != *otherNNeighsIter )
			return false;

		for ( const auto& thisToOuterDiff : *thisDiffIter )
			if( !MatchAny( thisToOuterDiff, *otherDiffIter ) )
				return false;
	}
	return true;
}

//---------------------------------------------------------------------------//

InterfaceCellIdentifierImp::InterfaceCellIdentifierImp (
	const PBCDifference* diffFunc
) noexcept :
	_DiffCorrection( diffFunc )
{}


pair<vector<Cell>, vector<Cell>>
InterfaceCellIdentifierImp::Identify (
	const vector<data3d>& innerPos,
	const vector<data3d>& outerPos
) const noexcept
{
	if ( innerPos.empty() ) return {vector<Cell>(), vector<Cell>()};

	auto inner2Outer = ParsingNeighboringStatus(innerPos, outerPos);
	cout<<"cellsInn2Out\n";
	auto inner2Inner = ParsingNeighboringStatus(innerPos, innerPos);
	cout<<"cellsInn2Inn\n";
/*int index=0;
for ( const auto& inn : inner2Inner )
{
cout<<index++<<": ";
for(unsigned i=0; i<inn.id.size(); ++i)
cout<<<"inn id:"<inn.id[i]<<" ";

cout<<endl;
cout<<"inn id size:"<<inn.id.size()<<endl;
}*/
	auto innerCells = IdentifyInnerCells( inner2Outer, inner2Inner );
	cout<<"auto innerCells\n";
	ChainingCells( innerCells, innerPos, innerCells.front(), innerPos );
	cout<<"ChainingCells\n";
	CheckInnerCellChaining( innerCells, innerPos, inner2Outer, inner2Inner );
	cout<<"Check\n";
cout<<"inner cells:"<<innerCells.size()<<"\n";
for(const auto& cell : innerCells)
{
cout<<cell.orientation<<": ";
for(auto id : cell.id)
cout<<id<<" ";
cout<<"("<<cell.position[0]<<" "<<cell.position[1]<<" "<<cell.position[2]<<")";
cout<<" - ["<<cell.x<<", "<<cell.y<<", "<<cell.z<<"]";
cout<<"\n";
}
cout<<endl;

	auto outerCells = IdentifyOuterCells( innerCells, innerPos, outerPos );
	ChainingCells( outerCells, outerPos, innerCells.front(), innerPos );
cout<<"outer cells:"<<outerCells.size()<<"\n";
for(const auto& cell : outerCells)
{
cout<<cell.orientation<<": ";
for(auto id : cell.id)
cout<<id<<" ";
cout<<"("<<cell.position[0]<<" "<<cell.position[1]<<" "<<cell.position[2]<<")";
cout<<" - ["<<cell.x<<", "<<cell.y<<", "<<cell.z<<"]";
cout<<"\n";
}
cout<<endl;
	return make_pair( move(innerCells), move(outerCells) );
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

vector<InterfaceCellIdentifierImp::Neighbors>
InterfaceCellIdentifierImp::ParsingNeighboringStatus (
	const vector<data3d>& from,
	const vector<data3d>& to
) const noexcept
{
	vector<Neighbors> neighs( from.size() );
	auto neighIter = neighs.begin();
	for (	auto fromIter = from.begin();
		fromIter != from.end();
		++fromIter, ++neighIter
	)
	{
		for (	auto toIter = to.begin();
			toIter != to.end();
			++toIter
		)
		{
			if ( fromIter == toIter ) continue;
			auto diff = Difference(*toIter, *fromIter);
			if ( OutOfRange(diff) ) continue;
			neighIter->id.push_back( distance(to.begin(), toIter) );
			neighIter->diff.push_back( diff );
		}
	}

	return neighs;
}


vector<Cell> 
InterfaceCellIdentifierImp::IdentifyInnerCells (
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner
) const noexcept
{
	unsigned nInners = inner2Outer.size();
	vector<bool> innerAtomIdentifyStatus( nInners, false );
	auto cells = IdentifyInnerCornerCells (
		inner2Outer, inner2Inner, innerAtomIdentifyStatus
	);
	std::cout<<"Corner\n";
	IdentifyInnerEdgeCells (
		inner2Outer, inner2Inner, innerAtomIdentifyStatus,
		cells
	);
	std::cout<<"Edge\n";
	IdentifyInnerFaceCells (
		inner2Outer, inner2Inner, innerAtomIdentifyStatus,
		cells
	);
	std::cout<<"Face\n";
	return cells;
}


vector<Cell>
InterfaceCellIdentifierImp::IdentifyOuterCells (
	const vector<Cell>& innerCells,
	const vector<data3d>& innerPos,
	const vector<data3d>& outerPos
) const noexcept
{
	vector<Cell> cells;
	for ( const auto& innerCell : innerCells )
	{
		auto outerCell = OuterCellOf( innerCell, innerPos, outerPos );
		cells.insert (
			cells.end(),
			make_move_iterator( outerCell.begin() ),
			make_move_iterator( outerCell.end() )
		);
	}
	return cells;
}

//---------------------------------------------------------------------------//

vector<Cell>
InterfaceCellIdentifierImp::IdentifyInnerCornerCells (
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner,
	vector<bool>& atomIdentifyStatus
) const noexcept
{
	vector<Cell> cells;
	for ( unsigned i=0; i<inner2Inner.size(); ++i )
	{
		if ( IsIdentified( i, atomIdentifyStatus ) ) continue;
		auto candidates = IdentifyInnerCellBasedOn (
			i, inner2Outer, inner2Inner, CORNER
		);
		std::cout<<candidates.size()<<"\n";
		if ( !candidates.empty() )
		{
			PrioritySort( candidates, atomIdentifyStatus );
			MarkCellIdentified( candidates[0], atomIdentifyStatus );
			cells.emplace_back( move(candidates[0]) );
		}
	}
	return cells;
}


void
InterfaceCellIdentifierImp::IdentifyInnerEdgeCells (
	const std::vector<Neighbors>& inner2Outer,
	const std::vector<Neighbors>& inner2Inner,
	std::vector<bool>& atomIdentifyStatus,
	std::vector<Cell>& cells
) const noexcept
{
	for ( unsigned i=0; i<inner2Inner.size(); ++i )
	{
		if ( IsIdentified( i, atomIdentifyStatus ) ) continue;
		auto candidates = IdentifyInnerCellBasedOn (
			i, inner2Outer, inner2Inner, EDGE
		);
		if ( !candidates.empty() )
		{
			DifferenceHint( i, cells, inner2Inner, candidates );
			PrioritySort( candidates, atomIdentifyStatus );
			MarkCellIdentified( candidates[0], atomIdentifyStatus );
			cells.emplace_back( move(candidates[0]) );
		}
	}
}


void
InterfaceCellIdentifierImp::IdentifyInnerFaceCells (
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner,
	std::vector<bool>& atomIdentifyStatus,
	std::vector<Cell>& cells
) const noexcept
{
	queue<pair<unsigned, CellOrientation>> seeds;
	AddSeeds( cells, seeds, inner2Inner );

	while (	!AllAtomsIdentified(seeds, atomIdentifyStatus) )
	{
		auto seed = seeds.front();
		seeds.pop();
		if ( IsIdentified( seed.first, atomIdentifyStatus ) ) continue;

		auto candidates = IdentifyInnerCellBasedOn (
			seed.first, inner2Outer, inner2Inner, FACE, seed.second
		);
		DifferenceHint( seed.first, cells, inner2Inner, candidates );
		PrioritySort ( candidates, atomIdentifyStatus );
		MarkCellIdentified( candidates[0], atomIdentifyStatus );
		AddSeeds( candidates[0], seeds, inner2Inner );
		cells.emplace_back( move(candidates[0]) );
	}
}

//---------------------------------------------//

vector<Cell>
InterfaceCellIdentifierImp::IdentifyInnerCellBasedOn (
	unsigned seed,
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner,
	CellType cellType,
	CellOrientation orientation
) const noexcept
{
	auto candidate = ProposeCandidates( seed, inner2Outer );
	PruneCandidates( candidate, inner2Outer, inner2Inner, orientation );
	std::cout<<candidate.size()<<"\n";
	auto endIter = candidate.end();
	for (	auto iter = candidate.begin();
		iter != endIter;
		++iter
	)
	{
		vector<Cell> thisCandidate;
		thisCandidate.reserve( candidate.capacity() );
		thisCandidate.emplace_back( *iter );
		while ( ExtendCandidatesCell( thisCandidate, inner2Inner ) )
			PruneCandidates( thisCandidate, inner2Outer, inner2Inner, orientation );
		candidate.insert (
			candidate.end(),
			make_move_iterator( thisCandidate.begin() ),
			make_move_iterator( thisCandidate.end() )
		);
	}
	std::cout<<"ITER DONE\n";
	candidate.erase( candidate.begin(), endIter );
	std::cout<<"erase all candidate\n";
	PruneCandidates( candidate, inner2Outer, inner2Inner, orientation, cellType );
	std::cout<<"Prune\n";
/*cout<<"total #candidates: "<<candidate.size()<<endl;
for(auto& c : candidate)
cout<<"@@"<<c.id[0]<<" "<<c.id[1]<<" "<<c.id[2]<<" "<<c.id[3]<<" "<<c.id[4]<<" "<<c.id[5]<<" "<<c.id[6]<<" "<<c.id[7]<<endl;
cout<<endl;*/
	return candidate;
}


vector<Cell>
InterfaceCellIdentifierImp::ProposeCandidates (
	unsigned orientation,
	const vector<Neighbors>& inner2Outer
) const noexcept
{
	auto N = NBasis();
	auto nInners = inner2Outer.size();

	vector<Cell> candidate;
	candidate.reserve( N*N*N );
	candidate.resize( N );
	for ( unsigned i=0; i<N; ++i )
	{
		candidate[i].id.assign( N, nInners+1 );
		candidate[i].id[i] = orientation;
	}
	return candidate;
}


bool
InterfaceCellIdentifierImp::ExtendCandidatesCell (
	vector<Cell>& candidate,
	const vector<Neighbors>& inner2Inner
) const noexcept
{
	auto nInner = inner2Inner.size();
	vector<Cell> proposeCandidate;
	for (	auto iter = candidate.begin();
		iter  != candidate.end();
		++iter
	)
	{
		auto mismatchIter = adjacent_find (
			iter->id.begin(), iter->id.end(),
			[nInner] (
				unsigned id1,
				unsigned id2
			) noexcept
			{
				return (id1>nInner) != (id2>nInner);
			}
		);
		if ( mismatchIter == iter->id.end() ) continue;

		auto mismatchIndex = distance( iter->id.begin(), mismatchIter );
		auto valueIndex = *mismatchIter>nInner ? mismatchIndex+1 : mismatchIndex;
		auto nullIndex = *mismatchIter>nInner ? mismatchIndex : mismatchIndex+1;
		if ( iter->id[valueIndex] == nInner || iter->id[nullIndex] == nInner ) continue;
		iter->id[nullIndex] = nInner;
		const auto& neighs = inner2Inner[ iter->id[valueIndex] ].id;
		const auto& neighDiff = inner2Inner[ iter->id[valueIndex] ].diff;
		for (	unsigned neighIndex=0; neighIndex<neighs.size(); ++neighIndex )
		{
			const auto& neighID = neighs[ neighIndex ];
			if ( !IsValidInnerRelation( valueIndex, nullIndex, neighDiff[neighIndex]) )
				continue;
			auto newCandidate( *iter );
			newCandidate.id[nullIndex] = neighID;
			proposeCandidate.emplace_back( move(newCandidate) );
		}
	}

	if ( proposeCandidate.empty() )
		return false;
	candidate.insert( candidate.end(), proposeCandidate.begin(), proposeCandidate.end() );
	return true;
}


void
InterfaceCellIdentifierImp::PruneCandidates (
	vector<Cell>& candidate,
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner,
	CellOrientation orientation,
	CellType cellType
) const noexcept
{

	auto newEnd = remove_if (
		candidate.begin(), candidate.end(),
		[this, &inner2Outer, &inner2Inner, cellType, orientation] (
			Cell& cell
		) noexcept
		{
			return !this->IsValidInnerCell (
				cell, inner2Outer, inner2Inner, orientation, cellType
			);
		}
	);
	candidate.erase( newEnd, candidate.end() );
}


InterfaceCellIdentifierImp::CellCondition
InterfaceCellIdentifierImp::ConvertToCellCondition (
	const Cell& cell,
	const vector<Neighbors>& inner2Outer
) const noexcept
{
	auto nInners = inner2Outer.size();
	auto nBasis = NBasis();

	CellCondition condition;
	auto& nNeigh2Other = condition.nNeigh2Other;
	auto& neighDiff = condition.neighDiff;
	nNeigh2Other.resize( nBasis );
	neighDiff.resize( nBasis );
	for ( unsigned i=0; i<nBasis; ++i )
	{
		const auto& id = cell.id[i];
		if ( id >= nInners ) continue;

		nNeigh2Other[i] = inner2Outer[id].id.size();
		neighDiff[i] = inner2Outer[id].diff;
	}

	return condition;
}


void
InterfaceCellIdentifierImp::PrioritySort (
	vector<Cell>& innerCells,
	const vector<bool>& atomIdentifyStatus
) const noexcept
{
	partial_sort (	innerCells.begin(), innerCells.begin()+1, innerCells.end(),
//	sort (	innerCells.begin(), innerCells.end(),
		[&atomIdentifyStatus] (
			const Cell& cellA,
			const Cell& cellB
		) noexcept
		{
			auto nInners = atomIdentifyStatus.size();

			auto foundID_A = 0;
			auto identified_A = false;
			for ( auto id : cellA.id )
			{
				if ( id < nInners )
					++foundID_A;
				identified_A |= atomIdentifyStatus[id];
			}
			auto foundID_B = 0;
			bool identified_B = false;
			for ( auto id : cellB.id )
			{
				if ( id < nInners )
					++foundID_B;
				identified_B |= atomIdentifyStatus[id];
			}
			if ( identified_A && identified_B )
				return Square(cellA.position) < Square(cellB.position);
			else if ( !identified_A && !identified_B )
				return foundID_A > foundID_B;
			else if ( identified_A )
				return false;
			else // identified_B
				return true;
		}
	);
}


void
InterfaceCellIdentifierImp::DifferenceHint (
	unsigned seed,
	const vector<Cell>& identifiedCells,
	const vector<Neighbors>& inner2Inner,
	vector<Cell>& candidateCells
) const noexcept
{
	const auto& refIDs = inner2Inner[seed].id;
	const auto& refDiffs = inner2Inner[seed].diff;
	data3d nearestDiff{};
	double nearestDiffSq = 10.;
	for ( const auto& cell : identifiedCells )
		for ( auto idIter=cell.id.begin(); idIter!=cell.id.end(); ++idIter )
		{
			auto iter = find(refIDs.begin(), refIDs.end(), *idIter);
			if ( iter != refIDs.end() )
			{
				auto dist = distance( refIDs.begin(), iter );
				auto diffIter = refDiffs.begin() + dist;
				auto diffSq = Square(*diffIter);
				if ( diffSq < nearestDiffSq )
				{
					nearestDiff = *diffIter;
					nearestDiffSq = diffSq;
				}
			}
		}

	for ( auto& cell : candidateCells )
	{
		double farestDiffSq = 0.;
		for ( auto id : cell.id )
		{
			auto refIDIter = find( refIDs.begin(), refIDs.end(), id );
			if ( refIDIter == refIDs.end() ) continue;

			auto refDiffIter = refDiffs.begin() + distance( refIDs.begin(), refIDIter );
			auto diff = Substract( *refDiffIter, nearestDiff );
			auto diffSq = Square( diff );
			if ( diffSq > farestDiffSq )
			{
				farestDiffSq = diffSq;
				cell.position = diff;
			}
		}
	}
}

//---------------------------------------------//

vector<Cell>
InterfaceCellIdentifierImp::OuterCellOf (
	const Cell& innerCell,
	const vector<data3d>& innerPos,
	const vector<data3d>& outerPos
) const noexcept
{
	vector<Cell> outerCells;
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, X, outerCells );
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, Y, outerCells );
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, Z, outerCells );
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, XY, outerCells );
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, XZ, outerCells );
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, YZ, outerCells );
	FindOuterCellAtoms ( innerCell, innerPos, outerPos, XYZ, outerCells );
	return outerCells;
}


void
InterfaceCellIdentifierImp::FindOuterCellAtoms (
	const Cell& innerCell,
	const vector<data3d>& innerPos,
	const vector<data3d>& outerPos,
	InterfaceCoordinate coord,
	vector<Cell>& outerCells
) const noexcept
{
	if ( !CellInterfaceAt(innerCell, coord) ) return;
	const auto& differenceToOuter = OuterCellDifference( coord, innerCell );

	auto nInners = innerPos.size();
	auto usedCellIter = find_if (
		innerCell.id.begin(), innerCell.id.end(),
		[nInners] ( unsigned id )
		{
			return id < nInners;
		}
	);
	auto usedIndex = distance( innerCell.id.begin(), usedCellIter );
	Cell outerCell = ReciprocalCell( innerCell, coord, outerPos.size() );
	for ( unsigned i=0; i<outerPos.size(); ++i )
	{
		auto difference = Difference( outerPos[i], innerPos[*usedCellIter] );
		unsigned j = MatchIndex(difference, differenceToOuter[usedIndex]);
		if ( j != differenceToOuter[usedIndex].size() )
			outerCell.id[j] = i;
	}
	outerCells.emplace_back( move(outerCell) );
}

//---------------------------------------------------------------------------//

void
InterfaceCellIdentifierImp::AddSeeds (
	const vector<Cell>& cells,
	queue<pair<unsigned, CellOrientation>>& seeds,
	const vector<Neighbors>& inner2Inner
) const noexcept
{
	for ( const auto& cell : cells )
		AddSeeds( cell, seeds, inner2Inner );
}


void
InterfaceCellIdentifierImp::AddSeeds (
	const Cell& cell,
	queue<pair<unsigned, CellOrientation>>& seeds,
	const vector<Neighbors>& inner2Inner
) const noexcept
{
	for ( auto id : cell.id )
		if ( id < inner2Inner.size() )
			for ( auto neigh : inner2Inner[id].id )
				seeds.push( make_pair(neigh, cell.orientation) );
}

//---------------------------------------------------------------------------//

bool
InterfaceCellIdentifierImp::AllAtomsIdentified (
	queue<pair<unsigned, CellOrientation>>& seeds,
	const vector<bool>& atomIdentifiedStatus
) const noexcept
{
	if ( !seeds.empty() ) return false;

	auto nonIdentifiedIter = find (
		atomIdentifiedStatus.begin(), atomIdentifiedStatus.end(),
		UNIDENTIFIED
	);
	if ( nonIdentifiedIter == atomIdentifiedStatus.end() )
		return true;

	seeds.push( make_pair (
		distance(atomIdentifiedStatus.begin(), nonIdentifiedIter),
		REGULAR // the first cell is assume as REGULAR arrangement
	) );
	return false;
}


void
InterfaceCellIdentifierImp::MarkCellIdentified (
	const Cell& cell,
	vector<bool>& atomIdentifiedStatus
) const noexcept
{
	for ( auto id : cell.id )
		if ( id < atomIdentifiedStatus.size() )
			atomIdentifiedStatus[id] = IDENTIFIED;
}


bool
InterfaceCellIdentifierImp::IsIdentified (
	unsigned seed,
	const vector<bool>& atomIdentifiedStatus
) const noexcept
{
	return atomIdentifiedStatus[seed];
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

bool
InterfaceCellIdentifierImp::OutOfRange (
	const data3d& diff
) const noexcept
{
	if ( fabs(diff[0])>1.00001 || fabs(diff[1])>1.00001 || fabs(diff[2])>1.00001 ) return true;
	auto diffSq = Square( diff );
	if ( diffSq > 1.00001 || diffSq < 0.00001 ) return true;
	return false;
}


data3d
InterfaceCellIdentifierImp::Difference (
	const data3d& to,
	const data3d& from
) const noexcept
{
	data3d diff = Substract( to, from );
	_DiffCorrection->Correction( diff.data() );
	return diff;
}

