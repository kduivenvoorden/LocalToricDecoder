/*
 * Automaton.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: michels
 */

#include "Automaton.h"
#include "AutomatonMemoryDataset.h"
#include "Location.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <list>

int Automaton::ipow(int base, int exp) {
    int result = 1;
    while (exp) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

Automaton::Automaton(Point arrayPosition, int L, int Q, int U, int b,
		int hierarchyDepth, double p, int fn, int fc, bool movie) :
		hierarchyDepth(hierarchyDepth), printMovie(movie){
	//this->fc = 3.0/5.0;
	//this->fn = 2.0/5.0;
	this->fn = fn; 
	this->fc = fc;
	for(int i = 0; i < 8; i++){
		this->syndromes[i] = 0;
	}
	for(int i = 1; i < hierarchyDepth + 1; i++){
		/*
		int colonyQ=0;
		if(i != hierarchyDepth){
			colonyQ = this->ipow(Q,i);
		} else {
			colonyQ = L;	
		}
		*/
		Point p = Point(arrayPosition.x % (this->ipow(Q,i)),
				arrayPosition.y % (this->ipow(Q,i)));
		this->memory.push_back(new AutomatonMemoryDataset(p,this->ipow(Q,i),
				this->ipow(U,i-1),this->ipow(b,i)));
		/*
		Point p = Point(arrayPosition.x % colonyQ,
				arrayPosition.y % colonyQ);
		this->memory.push_back(new AutomatonMemoryDataset(p,colonyQ,
				this->ipow(U,i),this->ipow(b,i)));
		*/
	}
}

Automaton::~Automaton() {
	delete[] this->neighbors;
	delete[] this->qubits;
	for(int i = 0; i < this->hierarchyDepth; i++){
		delete this->memory[i];
	}
}

void Automaton::setNeighbors(Automaton** neighbors) {
	this->neighbors = neighbors;
}

void Automaton::setQubits(Qubit** qubits) {
	this->qubits = qubits;
}

void Automaton::setSyndrome(int syndrome) {
	this->syndromes[Location::C] = syndrome;
}

void Automaton::printRule(int k, int depth){
	std::cout << depth <<" "<< this->memory[hierarchyDepth-1]->address.x<<" "<< this->memory[hierarchyDepth-1]->address.y<<" ";
	if (k==0)
		std::cout <<"1 0" <<std::endl;
	if (k==2)
		std::cout <<"0 1" <<std::endl;
	
	if (k==3)
		std::cout <<"-1 0" <<std::endl;
	if (k==1)
		std::cout <<"0 -1" <<std::endl;
}




int Automaton::localUpdateRules(Point address, int syndromes[9], Bit* bits[4],
		int Q) {
	
	
	
	// If no syndrome, do nothing
	if (syndromes[Location::C] == 0)
		return Location::C;
	
	// If at center, do nothing
	if( address.x == int(Q/2) && address.y == int(Q/2))
		return Location::C;
	
	// some ideas for improvement
	/*int tot = syndromes[Location::N]+syndromes[Location::E]+syndromes[Location::S]+syndromes[Location::W]+syndromes[Location::NE]+syndromes[Location::NW]+syndromes[Location::SE]+syndromes[Location::SW];
	if (tot==2){
		if (syndromes[Location::N] == 1 && (syndromes[Location::NE] == 1 || syndromes[Location::NW] == 1)){
			bits[Location::N]->flip();
			return Location::N;
		}
		if (syndromes[Location::E] == 1 && (syndromes[Location::NE] == 1 || syndromes[Location::SE] == 1)){
			bits[Location::E]->flip();
			return Location::E;
		}
		if (syndromes[Location::S] == 1 && (syndromes[Location::SE] == 1 || syndromes[Location::SW] == 1)){
			bits[Location::S]->flip();
			return Location::S;
		}
		if (syndromes[Location::W] == 1 && (syndromes[Location::NW] == 1 || syndromes[Location::SW] == 1)){
			bits[Location::W]->flip();
			return Location::W;
		}
	}*/
	
	// Nearest Neighbors
	std::list<int> NN_beAnnihalateBy;
	std::list<int> NN_canAnnihalate;

	
	if(address.x <= int(Q/2) && address.x > 0)
		NN_beAnnihalateBy.push_back(Location::W);
	if(address.x >= int(Q/2))
		NN_beAnnihalateBy.push_back(Location::E);
	if(address.y <= int(Q/2) && address.y > 0)
		NN_beAnnihalateBy.push_back(Location::S);
	if(address.y >= int(Q/2))
		NN_beAnnihalateBy.push_back(Location::N);
	
	
	
	//bourders
	if(address.x == 0)
		NN_canAnnihalate.push_back(Location::W);
	if(address.y == 0)
		NN_canAnnihalate.push_back(Location::S);
	

	
	if(address.x < int(Q/2))
		NN_canAnnihalate.push_back(Location::E);
	if(address.x > int(Q/2))
		NN_canAnnihalate.push_back(Location::W);
	if(address.y < int(Q/2))
		NN_canAnnihalate.push_back(Location::N);
	if(address.y > int(Q/2))
		NN_canAnnihalate.push_back(Location::S);
		
	
	
	for(int i : NN_beAnnihalateBy){
		if (syndromes[i] == 1)
			return Location::C;
	}
	for(int i : NN_canAnnihalate){
		if (syndromes[i] == 1){
			bits[i]->flip();
			return i;
		}
	}

	
	
	// Next Nearest Neighbors
	std::list<int> NNN_beAnnihalateBy;
	std::list<int> NNN_canAnnihalate;
	std::list<int> NNN_canAnnihalateAction;
	
	if(address.x <= int(Q/2) && address.y <= int(Q/2) && address.x > 0 && address.y > 0)
		NNN_beAnnihalateBy.push_back(Location::SW);
	if(address.x <= int(Q/2) && address.y >= int(Q/2) && address.x > 0 && address.y < Q-1)
		NNN_beAnnihalateBy.push_back(Location::NW);
	if(address.x >= int(Q/2) && address.y >= int(Q/2) && address.x <  Q-1 && address.y < Q-1)
		NNN_beAnnihalateBy.push_back(Location::NE);
	if(address.x >= int(Q/2) && address.y <= int(Q/2) && address.x < Q-1 && address.y > 0)
		NNN_beAnnihalateBy.push_back(Location::SE);

	//borders
	if(address.x == 0 && address.y <= int(Q/2)){
		NNN_canAnnihalate.push_back(Location::SW);
		NNN_canAnnihalateAction.push_back(Location::W);
	}
	if(address.x == 0 && address.y >= int(Q/2)){
		NNN_canAnnihalate.push_back(Location::NW);
		NNN_canAnnihalateAction.push_back(Location::W);
	}
	if(address.x <= int(Q/2) && address.y == 0){
		NNN_canAnnihalate.push_back(Location::SW);
		NNN_canAnnihalateAction.push_back(Location::S);
	}
	if(address.x <= int(Q/2) && address.y == Q-1){
		NNN_canAnnihalate.push_back(Location::NW);
		NNN_canAnnihalateAction.push_back(Location::N);
	}
	if(address.x >= int(Q/2) && address.y == 0){
		NNN_canAnnihalate.push_back(Location::SE);
		NNN_canAnnihalateAction.push_back(Location::S);
	}
	if(address.x >= int(Q/2) && address.y == Q-1){
		NNN_canAnnihalate.push_back(Location::NE);
		NNN_canAnnihalateAction.push_back(Location::N);
	}
	if(address.x == Q-1 && address.y <= int(Q/2)){
		NNN_canAnnihalate.push_back(Location::SE);
		NNN_canAnnihalateAction.push_back(Location::E);
	}
	if(address.x == Q-1 && address.y >= int(Q/2)){
		NNN_canAnnihalate.push_back(Location::NE);
		NNN_canAnnihalateAction.push_back(Location::E);
	}

	
	if(address.x <= int(Q/2) && address.y <= int(Q/2) ){
		NNN_canAnnihalate.push_back(Location::SE);
		NNN_canAnnihalateAction.push_back(Location::E);
		NNN_canAnnihalate.push_back(Location::NW);
		NNN_canAnnihalateAction.push_back(Location::N);
	}
	if(address.x <= int(Q/2) && address.y >= int(Q/2) ){
		NNN_canAnnihalate.push_back(Location::NE);
		NNN_canAnnihalateAction.push_back(Location::E);
		NNN_canAnnihalate.push_back(Location::SW);
		NNN_canAnnihalateAction.push_back(Location::S);
	}
	if(address.x >= int(Q/2) && address.y >= int(Q/2) ){
		NNN_canAnnihalate.push_back(Location::SE);
		NNN_canAnnihalateAction.push_back(Location::S);
		NNN_canAnnihalate.push_back(Location::NW);
		NNN_canAnnihalateAction.push_back(Location::W);
	}
	if(address.x >= int(Q/2) && address.y <= int(Q/2) ){
		NNN_canAnnihalate.push_back(Location::NE);
		NNN_canAnnihalateAction.push_back(Location::N);
		NNN_canAnnihalate.push_back(Location::SW);
		NNN_canAnnihalateAction.push_back(Location::W);
	}
	
	for(int i : NNN_beAnnihalateBy){
		if (syndromes[i] == 1)
			return Location::C;
	}
	while(NNN_canAnnihalate.size()>0){
		if (syndromes[NNN_canAnnihalate.front()] == 1){
			bits[NNN_canAnnihalateAction.front()]->flip();
			return NNN_canAnnihalateAction.front();
		}
		NNN_canAnnihalate.pop_front();
		NNN_canAnnihalateAction.pop_front();
	}
	
	
	
	
	// Remaining Action
	int move = Location::C;
	if( abs(address.x - int(Q/2)) > abs(address.y - int(Q/2))){
		if (address.x > int(Q/2))
			move = Location::W;
		if (address.x < int(Q/2))	
			move = Location::E;
	} else {	
		if (address.y > int(Q/2))
			move = Location::S;
		if (address.y < int(Q/2))	
			move = Location::N;
	}
	
	bits[move]->flip();
	return move;
	
}

void Automaton::localDependencies(){
	
	//Copies all information of neighbouring CA
	for(int i = 0; i < 8; i++){
		this->syndromes[i] = this->neighbors[i]->syndromes[Location::C];
	}
	
	// Increase age
	for(int i = 0; i < this->hierarchyDepth; i++)		
		this->memory[i]->age = (this->memory[i]->age + 1) % this->memory[i]->U;
	
	
	// TODO:: this copying is killing the speed do something
	//Updates all NEW- fields
	//newCountSignal
	//newisFlipedSignal
	//newFlipSignal
	for(int i = 0; i < this->hierarchyDepth - 1; i++){
		for(int j = 0; j < 8; j++){
			// move countSignal
			this->memory[i]->newCountSignal[Location::oppositeDirections[j]] =
					this->neighbors[j]->memory[i]->countSignal[Location::oppositeDirections[j]];
					
			// move flipSignal
			this->memory[i]->newIsFlipedSignal[Location::oppositeDirections[j]] =
				this->neighbors[j]->memory[i]->isFlipedSignal[Location::oppositeDirections[j]];

			// move isflipedSignal
			this->memory[i]->newFlipSignal[Location::oppositeDirections[j]]->bit =
				this->neighbors[j]->memory[i]->flipSignal[Location::oppositeDirections[j]]->bit;	
		}
		//At non colony centers, these will later be copied to the orginal fields, (non-new) fields. 
		//At     colony centers alternative calculations will be done to calculate (non-new) fields.
		//(at all times)
		//isFliped signal, from incoming newFlipsignal 
		//(end of workday)
		//CountSignal, at level 0 from syndrome, at level i, from syndromecount at level i-1
		//flipSignal, by local update rules
			
	}
}

void Automaton::executeFlips(){
	bool quickFlip = true;
	if(!quickFlip){
		for(int i = 0; i < this->hierarchyDepth-1; i++){
			if(this->memory[i+1]->age == 0){
				for(int j = 0; j < 4; j++){
					if(this->memory[i]->flipSignal[j]->bit == 1){
						this->qubits[j]->flip();
						this->memory[i]->flipSignal[j]->bit = 0;
						
					}
				}
			}
		}
	}
}

void Automaton::caUpdates(){
	// hard code the devision or non-devision strategy
	bool devision = false;
	bool quickFlip = true;
	int countBin = 0;
	int invert;
	Automaton* temp;
	std::list<int> binRange;
	if(devision){
		for(int i = 1; i <= this->memory[0]->b; i++)
			binRange.push_back(i);
	} else{
		binRange.push_back(0);
	}
	
	
	
	std::vector<int> centerSyndromes(9,0);
	memcpy(&centerSyndromes[0],&this->syndromes[0],9*sizeof(int));
	
	// THE CELLULAR AUTOMATON
	
	// On the first level it is a straightforward update rule. It flips bits depending on syndromes
	int k = this->localUpdateRules(this->memory[0]->address, &centerSyndromes[0], (Bit**)this->qubits, this->memory[0]->Q);
	if (k!=Location::C && printMovie)
		this->printRule(k,0);
	
	
	//On higher levels it is slightly more involved since one should first estimate the syndrome in the relavent representatives
	for(int i = 0; i < this->hierarchyDepth - 1; i++){
		
		// Colony center, representative, to ensure correct broadcasting, isFlipedSignal countSignal and flipSignal are calculated here .... 
		if(this->memory[i]->address.x == this->memory[i]->Q/2
				&& this->memory[i]->address.y == this->memory[i]->Q/2){
			
			// During workday
			if(this->memory[i]->age == 0){

				// countBin is used to seperate the countsignal into different bins.
				if (devision)
					countBin = (this->memory[i+1]->age/this->memory[i]->U)/this->memory[0]->b +1;
				
			
				// Update count and countSignal. 
				//For level 0 (at end of workday[-1], i.e. at all times)
				//count       depends on syndrome and newCountSingal
				//countSignal depends on syndrome 
				//For level i (at end of workday[i-1])
				//count       depends on count(i-1) and newCountSignal(i+1)
				//countSignal depends on count(i-1) 
				
				
				this->memory[i]->count[Location::C][countBin] += centerSyndromes.at(Location::C);			
				for(int j = 0; j < 8; j++){
					this->memory[i]->count[Location::getOppositeDirection(j)][countBin] += this->memory[i]->newCountSignal[j];
					this->memory[i]->countSignal[j] = centerSyndromes.at(Location::C);
				}
				
				// Note that centerSyndromes has not been reset and the value of the previous itteration is used (i.e. depending on count i)
				
				// Printing count signals at end of workday
				if(this->memory[i]->age == 0){	
					if (printMovie) {
						if (this->memory[i]->count[Location::C][0] > 0) {
							for( int l =0; l<9; l++)
								std::cout << -i-1 << " "<< this->memory[hierarchyDepth-1]->address.x<<" "<< this->memory[hierarchyDepth-1]->address.y << " "<< this->memory[i]->count[l][0] << " " << l<< std::endl;
						}
					}
				}

			
			}// During workday
			// End of workday
			if(this->memory[i+1]->age == 0){
				if (!devision) {
					invert = this->memory[1]->U;
				}else{
					invert = this->memory[0]->b;
				}
				
				
				// Inverting count depending on isFlippedSignal (i.e. compensate current conclusion for conclusions taken in previous end of workday)
				if(!quickFlip){
					for(int j = 0; j < 8; j++){
						if(this->memory[i]->newIsFlipedSignal[Location::getOppositeDirection(j)]==1){
							for(int l : binRange) 
								this->memory[i]->count[j][l] = invert - this->memory[i]->count[j][l];
						}
					}
					if (memory[i]->isFlipedSignal[Location::C]==1){
						for(int l : binRange){ 
							this->memory[i]->count[Location::C][l] = invert - this->memory[i]->count[Location::C][l];
							
						}
					}
				}

				if(devision){
					for(int j = 0; j < 9; j++){
						for(int l : binRange) {
							if(j == Location::C){
								if(this->memory[i]->count[j][l] >= this->fc)
									this->memory[i]->count[j][0] ++;
							} else{
								if(this->memory[i]->count[j][l] >= this->fn)
									this->memory[i]->count[j][0] ++;
							}
						}	
					}
				}
				
				
				
				
				// Estimates of syndromes of the colony centers
				// These will be used for local update rules as well as sending to higher levels.
				centerSyndromes.assign(9,0);

				//Calculate centerSyndromes (using fc and fn)
				for(int j = 0; j < 9; j++){	
					if(j == Location::C){
						if(this->memory[i]->count[j][0] >= this->fc)
							centerSyndromes.at(j) = 1;
					} else {
						if(this->memory[i]->count[j][0] >= this->fn)
							centerSyndromes.at(j) = 1;
					}
					
					// Reset count
					this->memory[i]->count[j].assign(this->memory[0]->b+1,0);
				}
				
				// determine the address of the center cell relative to the center cell of the higher hierarchy.
				// so to say factor the colonies out of the supercolonies (supercolonies/colonies) and take as
				// the representative for a colony the corresponding center cell
				int Q = this->memory[i]->Q;
				Point superAddress = this->memory[i+1]->address;
				Point reducedAddress = Point(superAddress.x/Q, superAddress.y/Q);
				
				//Rest data
				for(int j = 0; j < 9; j++)
					this->memory[i]->newFlipSignal[j]->bit = 0;
				this->memory[i]->isFlipedSignal.assign(9,0);
				
				// send corresponding flipSignals based on the local update rules
				int k = this->localUpdateRules(reducedAddress, &centerSyndromes[0], &this->memory[i]->flipSignal[0], this->memory[0]->Q);
				if (k!=Location::C){
					// this is used internaly to compensate in next workperiod (isFlipedSignal)
					// not that newFlipSignal[Location::C] is not reset in any other function
					this->memory[i]->newFlipSignal[Location::C]->bit = 1;
					if (printMovie) 
						this->printRule(k,i+1);
					
					
					if (quickFlip){
						this->qubits[k]->flip();
						temp = this->neighbors[k];
						for (int l = 1; l< this->memory[i]->Q;l++){
							temp->qubits[k]->flip();
							temp = temp->neighbors[k];
						}
					}
							
				}
				
				
			}// End of workday
			


			
			// The conclusion of a syndrome in a center should be adjust for possible planned flips. These should be re-calculated and re-sent at every time step
			//newFlipSignal is used since also incoming flipSignals of other colony centers are taken into account.
			int sum=0;
			for(int j = 0; j < 9; j++)
				sum += this->memory[i]->newFlipSignal[j]->bit;
			if(sum%2==1)
				this->memory[i]->isFlipedSignal.assign(9,1);
			
		} //colony centers
		// .... for all other cells, these are coppied straightforwardly from neighbouring information
		else {
			if(this->memory[i+1]->age == 0){
				//At the end of a workday should these be reset
				this->memory[i]->isFlipedSignal.assign(8,0);
				for(int j = 0; j < 8; j++)
					this->memory[i]->flipSignal[j]->bit = 0;
			} else {
				// otherwise these can be copied
				memcpy(&this->memory[i]->isFlipedSignal[0],&this->memory[i]->newIsFlipedSignal[0],8*sizeof(int));
				for(int j = 0; j < 8; j++)
				this->memory[i]->flipSignal[j]->bit = this->memory[i]->newFlipSignal[j]->bit;
			}
			
			// countSignals can always be copied
			memcpy(&this->memory[i]->countSignal[0],&this->memory[i]->newCountSignal[0],8*sizeof(int));
			if (printMovie){
				for(int j = 0; j < 8; j++){
					if (this->memory[i]->newCountSignal[j] > 0)
						std::cout << -i-10 << " "<< this->memory[hierarchyDepth-1]->address.x<<" "<< this->memory[hierarchyDepth-1]->address.y << " " << j<< std::endl;
				}
			}
		}
	}//hierarchy levels
}


/*int Automaton::localUpdateRules(Point address, int syndromes[9], Bit* bits[4],
		int Q) {
	// perform local update rules
	// w border
	if (address.x == 0) {
		if (syndromes[Location::C] == 0) {
			return Location::C;
		} else if (syndromes[Location::NW] == 1) {
			bits[Location::W]->flip();
			return Location::W;
		} else if (syndromes[Location::W] == 1) {
			bits[Location::W]->flip();
			return Location::W;
		} else if (syndromes[Location::SW] == 1) {
			bits[Location::W]->flip();
			return Location::W;
		}
	}
	// s border
	if (address.y == 0) {
		if (syndromes[Location::C] == 0) {
			return Location::C;
		}
		if (syndromes[Location::SW] == 1) {
			bits[Location::S]->flip();
			return Location::S;
		} else if (syndromes[Location::S] == 1) {
			bits[Location::S]->flip();
			return Location::S;
		} else if (syndromes[Location::SE] == 1) {
			bits[Location::S]->flip();
			return Location::S;
		}
	}
	if(address.x < int(Q/2)){
		// sw quadrant
		if (address.y < int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				bits[Location::N]->flip();
				return Location::N;
			} else if (syndromes[Location::E] == 1) {
				bits[Location::E]->flip();
				return Location::E;
			} else if (syndromes[Location::SW] == 1) {
				return Location::C;
			} else if (syndromes[Location::NW] == 1) {
				bits[Location::N]->flip();
				return Location::N;
			} else if (syndromes[Location::SE] == 1) {
				bits[Location::E]->flip();
				return Location::E;
			} else {
				//bits[Location::N]->flip();
				if(address.x - abs(address.y - Q) > 0){
					bits[Location::N]->flip();
					return Location::N;
				}
				else {
					bits[Location::E]->flip();
					return Location::E;
				}
			}
		}

		// w corridor
		if (address.y == int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				bits[Location::E]->flip();
				return Location::E;
			} else if (syndromes[Location::SW] == 1) {
				return Location::C;
			} else if (syndromes[Location::NW] == 1) {
				return Location::C;
			} else {
				bits[Location::E]->flip();
				return Location::E;
			}
		}

		// nw quadrant
		if (address.y > int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				bits[Location::E]->flip();
				return Location::E;
			} else if (syndromes[Location::S] == 1) {
				bits[Location::S]->flip();
				return Location::S;
			} else if (syndromes[Location::NW] == 1) {
				return Location::C;
			} else if (syndromes[Location::NE] == 1) {
				bits[Location::E]->flip();
				return Location::E;
			} else if (syndromes[Location::SW] == 1) {
				bits[Location::S]->flip();
				return Location::S;
			} else {
				if(address.x - address.y > 0){
					bits[Location::S]->flip();
					return Location::S;
				}else{
					bits[Location::E]->flip();
					return Location::E;
				}
			}
		}

	} else if (address.x > int(Q/2)) {
		// ne quadrant
		if (address.y > int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				bits[Location::S]->flip();
				return Location::S;
			} else if (syndromes[Location::W] == 1) {
				bits[Location::W]->flip();
				return Location::W;
			} else if (syndromes[Location::NE] == 1) {
				return Location::C;
			} else if (syndromes[Location::SE] == 1) {
				bits[Location::S]->flip();
				return Location::S;
			} else if (syndromes[Location::NW] == 1) {
				bits[Location::W]->flip();
				return Location::W;
			} else {
				if(abs(address.x-Q) - address.y > 0){
					bits[Location::S]->flip();
					return Location::S;
				}else{
					bits[Location::W]->flip();
					return Location::W;
				}
			}
		}
		// e corridor
		if (address.y == int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				bits[Location::W]->flip();
				return Location::W;
			} else if (syndromes[Location::NE] == 1) {
				return Location::C;
			} else if (syndromes[Location::SE] == 1) {
				return Location::C;
			} else {
				bits[Location::W]->flip();
				return Location::W;
			}
		}

		// se quadrant
		if (address.y < int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				bits[Location::W]->flip();
				return Location::W;
			} else if (syndromes[Location::N] == 1) {
				bits[Location::N]->flip();
				return Location::N;
			} else if (syndromes[Location::SE] == 1) {
				return Location::C;
			} else if (syndromes[Location::SW] == 1) {
				bits[Location::W]->flip();
				return Location::W;
			} else if (syndromes[Location::NE] == 1) {
				bits[Location::N]->flip();
				return Location::N;
			} else {
				if(abs(address.x-Q) - abs(address.y-Q) > 0){
					bits[Location::N]->flip();
					return Location::N;
				}else{
					bits[Location::W]->flip();
				//bits[Location::W]->flip();
					return Location::W;
				}
			}
		}
	} else {
		// n corridor
		if (address.y > int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				bits[Location::S]->flip();
				return Location::S;
			} else if (syndromes[Location::NW] == 1) {
				return Location::C;
			} else if (syndromes[Location::NE] == 1) {
				return Location::C;
			} else {
				bits[Location::S]->flip();
				return Location::S;
			}
		}


		// s corridor
		if (address.y < int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return Location::C;
			} else if (syndromes[Location::E] == 1) {
				return Location::C;
			} else if (syndromes[Location::S] == 1) {
				return Location::C;
			} else if (syndromes[Location::W] == 1) {
				return Location::C;
			} else if (syndromes[Location::N] == 1) {
				bits[Location::N]->flip();
				return Location::N;
			} else if (syndromes[Location::SE] == 1) {
				return Location::C;
			} else if (syndromes[Location::SW] == 1) {
				return Location::C;
			} else {
				bits[Location::N]->flip();
				return Location::N;
			}
		}
		return Location::C;
	}
}*/
