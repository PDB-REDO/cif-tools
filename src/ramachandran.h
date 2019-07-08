/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 19 juni, 2018
*/

#pragma once

float calculateRamachandranZScore(const std::string& aa, bool prePro, float phi, float psi);

enum RamachandranScore
{
	rsNotAllowed,
	rsAllowed,
	rsFavoured
};

RamachandranScore calculateRamachandranScore(const std::string& aa, bool prePro, float phi, float psi);

