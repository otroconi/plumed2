#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

int main() {

	std::vector<int> Bpair, myBpair;
	std::vector<double> bonddo, mybonddo;

	std::vector<int> Atrip, myAtrip;
	std::vector<double> thetao, mythetao;

	std::vector<int> Tquad, myTquad;
	std::vector<double> phio, myphio;

	std::vector<int> Ppairs, myPpairs;
	std::vector<double> pairsdo, mypairsdo;

	int check=1;

	size_t pos, p;
	string atoms, substr, at;
 
	size_t pos2, p2;
	string atoms2, substr2, at2;

	size_t pos3, p3;
	string atoms3, substr3, at3;

	size_t pos4, p4;
	string atoms4, substr4, at4;

	string line, line2, line3, line4;

	//File bonds.dat
	ifstream file ("bonds.dat");
	if (file.is_open())
	{
		while ( getline (file, line))
		{
			if (check%2 != 0) {
				pos = line.find("=");
				atoms = line.substr(pos+1);
				stringstream ss(atoms);
				while (ss.good()) {
					getline(ss, substr, ',');
					Bpair.push_back(std::stoi(substr));		
				}
	
//				cout << line << '\n';
			} else {
				pos = line.find("AT");
				p = 15;
				at = line.substr(pos+3,p);
				bonddo.push_back(std::stod(at,&p));		
//				cout << line << '\n';
			}
		++check;	
		}
		file.close();
	} 
	else cout << "Unable to open file bonds.dat" << endl;

/*
	for (int i=0; i<Bpair.size(); ++i) {
		cout << " " << Bpair[i] << endl;
	}


	for (int i=0; i<bonddo.size(); ++i) {
		cout << " " << bonddo[i] << endl;
	}


	ifstream inbp("Bpair.dat");

	if(!inbp) {
		std::cerr << "Cannot open the file: Bpair.dat" << std::endl;
		return false;
	}

	while ( getline(inbp, line)) {
		if (line.size() > 0)
			myBpair.push_back(std::stoi(line));

	}
	inbp.close();

	ifstream inbd("bonddo.dat");

	if(!inbd) {
		std::cerr << "Cannot open the file: bonddo.dat " << std::endl;
		return false;
	}

	while ( getline(inbd, line)) {
		if (line.size() > 0)
			mybonddo.push_back(std::stod(line));

	}
	inbd.close();

	
	for (int i=0; i<myBpair.size(); ++i) {
		cout << " " << myBpair[i] << endl;
	}

	for (int i=0; i<mybonddo.size(); ++i) {
		cout << " " << mybonddo[i] << endl;
	}
*/

	//File angles.dat
	ifstream file2 ("angles.dat");
	if (file2.is_open())
	{
		while ( getline (file2, line2))
		{
			if (check%2 != 0) {
				pos2 = line2.find("=");
				atoms2 = line2.substr(pos2+1);
				stringstream ss2(atoms2);
				while (ss2.good()) {
					getline(ss2, substr2, ',');
					Atrip.push_back(std::stoi(substr2));		
				}
	
//				cout << line2 << '\n';
			} else {
				pos2 = line2.find("AT");
				p2 = 15;
				at2 = line2.substr(pos2+3,p2);
				thetao.push_back(std::stod(at2,&p2));		
//				cout << line2 << '\n';
			}
		++check;	
		}
		file2.close();
	} 
	else cout << "Unable to open file angles.dat" << endl;
/*
	for (int i=0; i<Atrip.size(); ++i) {
		cout << " " << Atrip[i] << endl;
	}


	for (int i=0; i<thetao.size(); ++i) {
		cout << " " << thetao[i] << endl;
	}


	ifstream inat("Atrip.dat"); //file inat = input angle triplet

	if(!inat) {
		std::cerr << "Cannot open the file: Atrip.dat" << std::endl;
		return false;
	}

	while ( getline(inat, line2)) {
		if (line2.size() > 0)
			myAtrip.push_back(std::stoi(line2));

	}
	inat.close();

	ifstream inatheta("thetao.dat"); //file inad = input angle theta

	if(!inatheta) {
		std::cerr << "Cannot open the file: thetao.dat " << std::endl;
		return false;
	}

	while ( getline(inatheta, line2)) {
		if (line2.size() > 0)
			mythetao.push_back(std::stod(line2));

	}
	inatheta.close();

	
	for (int i=0; i<myAtrip.size(); ++i) {
		cout << " " << myAtrip[i] << endl;
	}

	for (int i=0; i<mythetao.size(); ++i) {
		cout << " " << mythetao[i] << endl;
	}
	
//	cout << " " << mythetao.size() << endl;

*/	
	//File torsions.dat
	ifstream file3 ("torsions.dat");
	if (file3.is_open())
	{
		while ( getline (file3, line3))
		{
			if (check%2 != 0) {
				pos3 = line3.find("=");
				atoms3 = line3.substr(pos3+1);
				stringstream ss3(atoms3);
				while (ss3.good()) {
					getline(ss3, substr3, ',');
					Tquad.push_back(std::stoi(substr3));		
				}
	
//				cout << line3 << '\n';
			} else {
				pos3 = line3.find("AT");
				p3 = 15;
				at3 = line3.substr(pos3+3,p3);
				phio.push_back(std::stod(at3,&p3));		
//				cout << line3 << '\n';
			}
		++check;	
		}
		file3.close();
	} 
	else cout << "Unable to open file torsions.dat" << endl;

/*
	for (int i=0; i<Tquad.size(); ++i) {
		cout << " " << Tquad[i] << endl;
	}


	for (int i=0; i<phio.size(); ++i) {
		cout << " " << phio[i] << endl;
	}


	ifstream intq("Tquad.dat"); //file intq = input torsion quadruple

	if(!intq) {
		std::cerr << "Cannot open the file: Tquad.dat" << std::endl;
		return false;
	}

	while ( getline(intq, line3)) {
		if (line3.size() > 0)
			myTquad.push_back(std::stoi(line3));

	}
	intq.close();

	ifstream intphi("phio.dat"); //file intphi = input torsion phi

	if(!intphi) {
		std::cerr << "Cannot open the file: phio.dat " << std::endl;
		return false;
	}

	while ( getline(intphi, line3)) {
		if (line3.size() > 0)
			myphio.push_back(std::stod(line3));

	}
	intphi.close();

	
	for (int i=0; i<myTquad.size(); ++i) {
		cout << " " << myTquad[i] << endl;
	}

	for (int i=0; i<myphio.size(); ++i) {
		cout << " " << myphio[i] << endl;
	}
	
	cout << " " << myphio.size() << endl;

*/

	//File pairs.dat
	ifstream file4 ("pairs.dat");
	if (file4.is_open())
	{
		while ( getline (file4, line4))
		{
			if (check%2 != 0) {
				pos4 = line4.find("=");
				atoms4 = line4.substr(pos4+1);
				stringstream ss4(atoms4);
				while (ss4.good()) {
					getline(ss4, substr4, ',');
					Ppairs.push_back(std::stoi(substr4));		
				}
	
//				cout << line4 << '\n';
			} else {
				pos4 = line4.find("FUNC");
				p4 = 15;
				at4 = line4.substr(pos4+17,p4);
				pairsdo.push_back(std::stod(at4,&p4));		
//				cout << line4 << '\n';
			}
		++check;	
		}
		file4.close();
	} 
	else cout << "Unable to open file pairs.dat" << endl;

/*
	for (int i=0; i<Ppairs.size(); ++i) {
		cout << " " << Ppairs[i] << endl;
	}


	for (int i=0; i<pairsdo.size(); ++i) {
		cout << " " << pairsdo[i] << endl;
	}


	ifstream inpp("Ppairs.dat");

	if(!inpp) {
		std::cerr << "Cannot open the file: Ppairs.dat" << std::endl;
		return false;
	}

	while ( getline(inpp, line4)) {
		if (line4.size() > 0)
			myPpairs.push_back(std::stoi(line4));

	}
	inpp.close();

	ifstream inpd("pairsdo.dat");

	if(!inpd) {
		std::cerr << "Cannot open the file: pairsdo.dat " << std::endl;
		return false;
	}

	while ( getline(inpd, line4)) {
		if (line4.size() > 0)
			mypairsdo.push_back(std::stod(line4));

	}
	inpd.close();
	
	for (int i=0; i<myPpairs.size(); ++i) {
		cout << " " << myPpairs[i] << endl;
	}

	for (int i=0; i<mypairsdo.size(); ++i) {
		cout << " " << mypairsdo[i] << endl;
	}
*/
	return 0;

}
