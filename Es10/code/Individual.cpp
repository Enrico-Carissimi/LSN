#include "Individual.hpp"



Individual::Individual(int l, Random& random){
    size = l;
    score = 0.;

    // the first element is fixed to 1, generate the other l-1 randomly
    vector<int> s;
    for (int i = 0; i < size - 1; i++){s.push_back(i);} // ordered sequence of 0, ..., l-2
    int currentL = size - 1;

    while (!s.empty()){
        int j = random.Rannyu() * currentL; // index between 0 and [current size of s] - 1
        genes.push_back(s[j] + 2);          // add random element from s to genes (+2 since we need 2, ..., l)

        s.erase(s.begin() + j);             // erase used elements to avoid repetition
        currentL--;                         // size of s is reduced by 1
    }

    genes.insert(genes.begin(), 1);         // insert 1 as first element

    check();
}



Individual Individual::operator=(const Individual& other){
    genes = other.getGenes();
    size = other.getSize();
    score = other.getScore();

    return *this;
}



void Individual::swapPair(int i, int j){
    swap(genes[i], genes[j]);
    check();
}

void Individual::swapPair(Random& random){
    int i = random.Integer(1, size), j = random.Integer(1, size); // from 1 to l-1
    
    swapPair(i, j);

    check();
}



void Individual::swapNeighbours(Random& random){
    int i = random.Integer(1, size);
    int neighbour = (random.Rannyu() < 0.5) ? 1 : -1;

    int j = i + neighbour;
    if (i == size - 1) j = i - 1;
    if (i == 1) j = i + 1;

    swapPair(i, j);
}



void Individual::flipBlock(int start, int blockL){
    if (start + blockL > size){
        cout << "Illegal action: trying to access out of bounds element" << endl;
        exit(-1);
    }

    // if blockL is odd, blockL/2 is rounded to (blockL-1)/2, thus this cicle won't affect the central element of the block
    for (int i = 0; i < blockL / 2; i++){
        swapPair(start + i, start + blockL - 1 - i);
    }

    check();
}

void Individual::flipBlock(Random& random){
    int start = random.Integer(1, size);
    int l = random.Integer(1, size - start);

    flipBlock(start, l);
}



void Individual::swapBlock(int start, int blockL, int destination){
    if (start + blockL > size || destination + blockL >  size){
        cout << "Illegal action: trying to access out of bounds element" << endl;
        exit(-1);
    }
    if (destination < start + blockL && destination >= start){
        cout << "Illegal action: trying to move block inside itself (swap)" << endl;
        exit(-1);
    }

    for (int i = 0; i < blockL; i++){
        swapPair(start + i, destination + i);
    }

    check();
}

void Individual::swapBlock(Random& random){
    // start or destination cannot be the first element since it is fixed
    int start = random.Integer(1, size);

    int destination;
    do{
        destination = random.Integer(1, size);
    }
    while (destination == start);

    if (destination < start) swap(destination, start);

    int maxL = min(destination - start, size - destination); 
    int l = random.Integer(1, maxL);

    swapBlock(start, l, destination);
}



void Individual::moveBlock(int start, int blockL, int destination){
    if (start + blockL > size){
        cout << "Illegal action: trying to access out of bounds element" << endl;
        exit(-1);
    }
    if (destination > size){
        cout << "Illegal action: trying to move block out of bounds (move)" << endl;
        exit(-1);
    }
    if (destination <= start + blockL && destination >= start){
        cout << "Illegal action: trying to move block inside itself" << endl;
        exit(-1);
    }

    vector<int> buffer;

    for (int i = 0; i < blockL; i++){
        buffer.push_back(genes[start + i]);
    }

    genes.insert(genes.begin() + destination, buffer.begin(), buffer.end());
    
    if (start < destination) genes.erase(genes.begin() + start, genes.begin() + start + blockL);
    else genes.erase(genes.begin() + blockL + start, genes.begin() + start + 2 * blockL);

    check();
}

void Individual::moveBlock(Random& random){
    // start cannot be the first element since it is fixed
    int start = random.Integer(1, size);
    int l = random.Integer(1, size - start);
    
    // generate a destination outside the selected block (and is not the first element)
    int destination = random.Integer(1, size - l);
    if (destination >= start) destination += l + 1;

    moveBlock(start, l, destination);
}



int Individual::find(int n){
    if (n > size || n < 1){
        cout << "Illegal action: searching for " << n << " (size is " << size << ")" << endl;
        exit(-1);
    }
    bool found = false;
    int i = 0;
    while (!found){
        if (genes[i] == n) found = true;
        else i++;
    }
    return i;
}



void Individual::check(){
    if (genes[0] != 1){
        cout << "Illegal individual: first city is wrong (not 1)" << endl;
        cout << "\t";
        print();
        exit(-1);
    }

    if ((int)genes.size() != size){
        cout << "Illegal individual: inconsistent size" << endl;
        cout << "\t";
        print();
        exit(-1);
    }

    // if a number has not yet appeared, set the corresponding element of c to true
    // if a has already appeared, it means that it is duplicated
    // if a number is missing and the size is correct, another one must be duplicated
    // this also checks if some number is outside the range (1, size)
    vector<bool> c(size, false);
    for (int i = 0; i < size; i++){
        int j = genes[i] - 1; // genes start from 1

        if (j > size){
            cout << "Illegal individual: value too high (" << j + 1 << ")" << endl;
            cout << "\t";
            print();
            exit(-1);
        }
        else if (j < 0){
            cout << "Illegal individual: value too low (" << j + 1<< ")" << endl;
            cout << "\t";
            print();
            exit(-1);
        }

        if (!c[j]) c[j] = true;
        else{
            cout << "Illegal individual: missing/duplicated city: " << j + 1 << endl;
            cout << "\t";
            print();
            exit(-1);
        }
    }
}



void Individual::print(){
    cout << "(";
    for (int i = 0; i < size; i++){
        cout << genes[i] << ", ";
    }
    cout << "\b\b)" << endl;
}

void Individual::print(string fileName){
    ofstream out(fileName, ios::app);

    for (int i = 0; i < size; i++){
        out << setw(4) << genes[i];
    }
    out << "\t" << getScore() << endl;

    out.close();
}