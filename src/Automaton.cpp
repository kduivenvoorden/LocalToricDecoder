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
		int hierarchyDepth, double p, int fn, int fc) :
		hierarchyDepth(hierarchyDepth){
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
				this->ipow(U,i),this->ipow(b,i)));
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

void Automaton::localUpdateRules(Point address, int syndromes[9], Bit* bits[4],
		int Q) {
	// perform local update rules
	// w border
	if (address.x == 0) {
		if (syndromes[Location::C] == 0) {
			return;
		} else if (syndromes[Location::NW] == 1) {
			bits[Location::W]->flip();
			return;
		} else if (syndromes[Location::W] == 1) {
			bits[Location::W]->flip();
			return;
		} else if (syndromes[Location::SW] == 1) {
			bits[Location::W]->flip();
			return;
		}
	}
	// s border
	if (address.y == 0) {
		if (syndromes[Location::C] == 0) {
			return;
		}
		if (syndromes[Location::SW] == 1) {
			bits[Location::S]->flip();
			return;
		} else if (syndromes[Location::S] == 1) {
			bits[Location::S]->flip();
			return;
		} else if (syndromes[Location::SE] == 1) {
			bits[Location::S]->flip();
			return;
		}
	}
	if(address.x < int(Q/2)){
		// sw quadrant
		if (address.y < int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::S] == 1) {
				return;
			} else if (syndromes[Location::W] == 1) {
				return;
			} else if (syndromes[Location::N] == 1) {
				bits[Location::N]->flip();
				return;
			} else if (syndromes[Location::E] == 1) {
				bits[Location::E]->flip();
				return;
			} else if (syndromes[Location::SW] == 1) {
				return;
			} else if (syndromes[Location::NW] == 1) {
				bits[Location::N]->flip();
				return;
			} else if (syndromes[Location::SE] == 1) {
				bits[Location::E]->flip();
				return;
			} else {
				//bits[Location::N]->flip();
				if(address.x - abs(address.y - Q) > 0)
					bits[Location::N]->flip();
				else
					bits[Location::E]->flip();
				return;
			}
		}

		// w corridor
		if (address.y == int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::S] == 1) {
				return;
			} else if (syndromes[Location::W] == 1) {
				return;
			} else if (syndromes[Location::N] == 1) {
				return;
			} else if (syndromes[Location::E] == 1) {
				bits[Location::E]->flip();
				return;
			} else if (syndromes[Location::SW] == 1) {
				return;
			} else if (syndromes[Location::NW] == 1) {
				return;
			} else {
				bits[Location::E]->flip();
				return;
			}
		}

		// nw quadrant
		if (address.y > int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::W] == 1) {
				return;
			} else if (syndromes[Location::N] == 1) {
				return;
			} else if (syndromes[Location::E] == 1) {
				bits[Location::E]->flip();
				return;
			} else if (syndromes[Location::S] == 1) {
				bits[Location::S]->flip();
				return;
			} else if (syndromes[Location::NW] == 1) {
				return;
			} else if (syndromes[Location::NE] == 1) {
				bits[Location::E]->flip();
				return;
			} else if (syndromes[Location::SW] == 1) {
				bits[Location::S]->flip();
				return;
			} else {
				if(address.x - address.y > 0)
					bits[Location::S]->flip();
				else
					bits[Location::E]->flip();
				//bits[Location::E]->flip();
				return;
			}
		}

	} else if (address.x > int(Q/2)) {
		// ne quadrant
		if (address.y > int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::N] == 1) {
				return;
			} else if (syndromes[Location::E] == 1) {
				return;
			} else if (syndromes[Location::S] == 1) {
				bits[Location::S]->flip();
				return;
			} else if (syndromes[Location::W] == 1) {
				bits[Location::W]->flip();
				return;
			} else if (syndromes[Location::NE] == 1) {
				return;
			} else if (syndromes[Location::SE] == 1) {
				bits[Location::S]->flip();
				return;
			} else if (syndromes[Location::NW] == 1) {
				bits[Location::W]->flip();
				return;
			} else {
				if(abs(address.x-Q) - address.y > 0)
					bits[Location::S]->flip();
				else
					bits[Location::W]->flip();
				//bits[Location::S]->flip();
				return;
			}
		}
		// e corridor
		if (address.y == int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::N] == 1) {
				return;
			} else if (syndromes[Location::E] == 1) {
				return;
			} else if (syndromes[Location::S] == 1) {
				return;
			} else if (syndromes[Location::W] == 1) {
				bits[Location::W]->flip();
				return;
			} else if (syndromes[Location::NE] == 1) {
				return;
			} else if (syndromes[Location::SE] == 1) {
				return;
			} else {
				bits[Location::W]->flip();
				return;
			}
		}

		// se quadrant
		if (address.y < int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::E] == 1) {
				return;
			} else if (syndromes[Location::S] == 1) {
				return;
			} else if (syndromes[Location::W] == 1) {
				bits[Location::W]->flip();
				return;
			} else if (syndromes[Location::N] == 1) {
				bits[Location::N]->flip();
				return;
			} else if (syndromes[Location::SE] == 1) {
				return;
			} else if (syndromes[Location::SW] == 1) {
				bits[Location::W]->flip();
				return;
			} else if (syndromes[Location::NE] == 1) {
				bits[Location::N]->flip();
				return;
			} else {
				if(abs(address.x-Q) - abs(address.y-Q) > 0)
					bits[Location::N]->flip();
				else
					bits[Location::W]->flip();
				//bits[Location::W]->flip();
				return;
			}
		}
	} else {
		// n corridor
		if (address.y > int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::W] == 1) {
				return;
			} else if (syndromes[Location::N] == 1) {
				return;
			} else if (syndromes[Location::E] == 1) {
				return;
			} else if (syndromes[Location::S] == 1) {
				bits[Location::S]->flip();
				return;
			} else if (syndromes[Location::NW] == 1) {
				return;
			} else if (syndromes[Location::NE] == 1) {
				return;
			} else {
				bits[Location::S]->flip();
				return;
			}
		}


		// s corridor
		if (address.y < int(Q / 2)) {
			if (syndromes[Location::C] == 0) {
				return;
			} else if (syndromes[Location::E] == 1) {
				return;
			} else if (syndromes[Location::S] == 1) {
				return;
			} else if (syndromes[Location::W] == 1) {
				return;
			} else if (syndromes[Location::N] == 1) {
				bits[Location::N]->flip();
				return;
			} else if (syndromes[Location::SE] == 1) {
				return;
			} else if (syndromes[Location::SW] == 1) {
				return;
			} else {
				bits[Location::N]->flip();
				return;
			}
		}
	}
}

void Automaton::localUpdates(){
	for(int i = 0; i < 8; i++){
		this->syndromes[i] = this->neighbors[i]->syndromes[Location::C];
	}
	this->localUpdateRules(this->memory[0]->address, this->syndromes, (Bit**)this->qubits, this->memory[0]->Q);

	for(int i = 0; i < this->hierarchyDepth; i++){
		// increment age
		this->memory[i]->age = (this->memory[i]->age + 1) % this->memory[i]->U;

		// TODO:: this copying is killing the speed do something
		// move countSignal
		for(int j = 0; j < 8; j++){
			this->memory[i]->newCountSignal[Location::oppositeDirections[j]] =
					this->neighbors[j]->memory[i]->countSignal[Location::oppositeDirections[j]];
		}
		// move flipSignal
		for(int j = 0; j < 8; j++){
			this->memory[i]->newFlipSignal[Location::oppositeDirections[j]]->bit =
				this->neighbors[j]->memory[i]->flipSignal[Location::oppositeDirections[j]]->bit;
		}
		// calculate new count of syndromes at the center of a colony and finalize movements by copying from
		// temporary storage into real storage
		if(this->memory[i]->address.x == this->memory[i]->Q/2
				&& this->memory[i]->address.y == this->memory[i]->Q/2){
			for(int j = 0; j < 8; j++)
				this->memory[i]->countSignal[j] = this->syndromes[Location::C];

			for(int j = 0; j < 9; j++){
				if(j == Location::C){
					if(this->syndromes[j] != 0)
						this->memory[i]->count[j][0]++;
				} else {
					if (this->memory[i]->newCountSignal[j] != 0){
						this->memory[i]->count[Location::getOppositeDirection(j)][0]++;
					}
				}
			}

			// uncomment the next section if you want to use the b subdivision for the count field
			/*
			if (this->memory[i]->age % this->memory[i]->b == 0){
				for(int j = 0; j < 8; j++){
					//if(this->memory[i]->count[j][0] >= this->fn * this->memory[i]->b){
					if(this->memory[i]->count[j][0] >= this->fn){
						this->memory[i]->count[j][0] = 0;
						this->memory[i]->count[j][1]++;
					}
				}
				//if(this->memory[i]->count[Location::C][0] >= this->fc * this->memory[i]->b){
				if(this->memory[i]->count[Location::C][0] >= this->fc){
					this->memory[i]->count[Location::C][0] = 0;
					this->memory[i]->count[Location::C][1]++;
				}
			}
			*/
			
		} else {
			memcpy(&this->memory[i]->countSignal[0],&this->memory[i]->newCountSignal[0],8*sizeof(int));
			for(int j = 0; j < 8; j++){
				this->memory[i]->flipSignal[j]->bit = this->memory[i]->newFlipSignal[j]->bit;
			}
		}
	}
}

void Automaton::executeFlips(){
	for(int i = 0; i < this->hierarchyDepth-1; i++){
		if(this->memory[i]->age % this->memory[i]->U == 0){
			for(int j = 0; j < 8; j++){
				if(this->memory[i]->flipSignal[j]->bit == 1){
					this->qubits[j]->flip();
					// we need to save the previous flipSignal in order to determine later (when the flipSignal field
					// has been reset) whether to set a new flipSignal. it's convenient to use the newFlipSignal
					// field for that purpose
					this->memory[i]->newFlipSignal[j]->bit = this->memory[i]->flipSignal[j]->bit;
					// reset flip signal
					this->memory[i]->flipSignal[j]->bit = 0;
				}
			}
		}
	}
}

void Automaton::signalUpdates(){
	// now we need to perform the 'local' updates on the hierarchy level, so to say pair out syndromes at colony,
	// supercolony,... centers
	// we skip the highestHierarchy because it includes the whole lattice
	for(int i = 0; i < this->hierarchyDepth - 1; i++){
		// updates are only necessary when the workday in a hierarchy level is over
		if(this->memory[i]->age % this->memory[i]->U == 0){
			if(this->memory[i]->address.x == this->memory[i]->Q/2
					&& this->memory[i]->address.y == this->memory[i]->Q/2){
				std::vector<int> centerSyndromes(9,0);

				for(int j = 0; j < 9; j++){
					// uncomment the next section if you want to use the b subdivision 
					/*
					if(j == Location::C){
						//if(this->memory[i]->count[j][1] >= this->fc * this->memory[i]->b)
						if(this->memory[i]->count[j][1] >= this->fc)
							centerSyndromes[j] = 1;
					} else {
						//if(this->memory[i]->count[j][1] >= this->fn * this->memory[i]->b)
						if(this->memory[i]->count[j][1] >= this->fn)
							centerSyndromes[j] = 1;
					}
					*/
					
					// comment the next section if you want to use the b subdivision
					if(i == 0){
						if(j == Location::C){
							if(this->memory[i]->count[j][0] >= this->fc)
								centerSyndromes[j] = 1;
						} else {
							if(this->memory[i]->count[j][0] >= this->fn)
								centerSyndromes[j] = 1;
						}
					} else {
						if(j == Location::C){
							if(this->memory[i]->count[j][0]/this->memory[i-1]->U >= this->fc)
								centerSyndromes[j] = 1;
						} else {
							if(this->memory[i]->count[j][0]/this->memory[i-1]->U >= this->fn)
								centerSyndromes[j] = 1;
						}
					}
					
					// the work of count is done for this day, reset it
					// TODO: maybe only reset when syndromes is set
					this->memory[i]->count[j][0] = 0;
					this->memory[i]->count[j][1] = 0;
				}
				// determine the address of the center cell relative to the center cell of the higher hierarchy.
				// so to say factor the colonies out of the supercolonies (supercolonies/colonies) and take as
				// the representative for a colony the corresponding center cell
				int Q = this->memory[i]->Q;
				Point superAddress = this->memory[i+1]->address;
				Point reducedAddress = Point(superAddress.x/Q, superAddress.y/Q);
				// determine if flip signal has already happened in the 'executeFlip' step and therefore do not
				// send a flipSignal again

				int sum = 0;
				for(int j = 0; j < 8; j++){
					centerSyndromes[j] = centerSyndromes[j] ^ this->memory[i]->newFlipSignal[j]->bit;
					sum += this->memory[i]->newFlipSignal[j]->bit;
				}
				sum = sum % 2;
				centerSyndromes[Location::C] = centerSyndromes[Location::C] ^ (sum);

				//int reducedQ = this->memory[i+1]->Q/this->memory[i]->Q;
				// send corresponding flipSignals based on the local update rules
				this->localUpdateRules(reducedAddress, &centerSyndromes[0], &this->memory[i]->flipSignal[0], this->memory[0]->Q);
				//this->localUpdateRules(reducedAddress, &centerSyndromes[0], &this->memory[i]->flipSignal[0], reducedQ);
			}
		}
	}
}
