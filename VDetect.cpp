#include "vdetect.h"
VDetect::VDetect(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    //checks if the size is prime
    if (isPrime(size)){
        //checks if size is less than minPrime
        if (size < MINPRIME){
            size = MINPRIME;
        }
        //checks if size is greater than maxPrime
        else if (size > MAXPRIME){
            size = MAXPRIME;
        }
    }
    //checks if the size it not prime and calls findNextPrime
    else if (!isPrime(size)){
        findNextPrime(size);
    }

    //initializes member variables
    m_currentCap = size;
    m_currentSize = 0;
    m_hash = hash;
    m_newPolicy = probing;
    m_currNumDeleted = 0;
    m_currProbing = probing;

    m_currentTable = new Virus[m_currentCap];

    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    m_oldProbing = probing;

}

VDetect::~VDetect(){
    //deletes the current table and old table
    delete [] m_currentTable;
    delete [] m_oldTable;

    m_currentCap = 0;
    m_currentSize = 0;
    m_currNumDeleted = 0;

    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;

}

void VDetect::changeProbPolicy(prob_t policy){
    //changes the probing policy
    m_newPolicy = policy;

}


bool VDetect::insert(Virus virus){

    Virus tempVirus = getVirus(virus.m_key, virus.m_id);
    //check if the virus already exists
    if (tempVirus.m_key != "" and tempVirus.m_id != -1){
        return false;
    }

    //calls the insert helper
    insertHelper(virus);

    //checks if the load factor is greater than 0.5
    float loadFactor = lambda();
    if ((loadFactor > 0.5) or (m_oldTable != nullptr)){
        rehash();
    }

    return true;
}

void VDetect::insertHelper(Virus virus){
    //checks if the virus is valid
    if (virus.m_id < MINID || virus.m_id > MAXID){
        return;
    }

    //creates i and sets index to the hash mod capacity
    int i = 0;
    int index = m_hash(virus.m_key) % m_currentCap;

    //iterates while the table at index is not empty or deleted
    while (!(m_currentTable[index] == EMPTY or (m_currentTable[index] == DELETED))){
        if (m_currProbing == NONE){
            return;
        }

            //calculates the index for quadratic
        else if (m_currProbing == QUADRATIC){
            index = ((m_hash(virus.m_key) % m_currentCap) + (i * i)) % m_currentCap;
        }

            //calculates the index for doubleHash
        else if (m_currProbing == DOUBLEHASH){
            index = ((m_hash(virus.m_key) % m_currentCap) + i * (11 - (m_hash(virus.m_key) % 11))) % m_currentCap;
        }
        i++;
    }

    //sets the table at index to virus and increases the size
    m_currentTable[index] = virus;
    m_currentSize++;

}

void VDetect::rehash(){

    if (m_oldTable == nullptr){
        //Save the old table information
        m_oldTable = m_currentTable;
        m_oldCap = m_currentCap;
        m_oldSize = m_currentSize;
        m_oldNumDeleted = m_currNumDeleted;
        m_oldProbing = m_currProbing;
        m_currProbing = m_newPolicy;

        // Find the next prime number for the new table size
        m_currentCap= findNextPrime(((m_currentSize - m_oldNumDeleted)*4));

        //create new table
        m_currentTable = new Virus[m_currentCap];
        for (int i = 0; i < m_currentCap; i++){
            m_currentTable[i] = EMPTY;
        }
        m_currentSize = 0;
        m_currNumDeleted = 0;
        m_currProbing = m_oldProbing;
    }

    //performs transfer
    int numTransferred = 0;
    for (int i = 0; i < m_oldCap && numTransferred < m_oldSize * 0.25; i++){
        if (!(m_oldTable[i] == EMPTY) && !(m_oldTable[i] == DELETED)){
            insertHelper(m_oldTable[i]);
            numTransferred++;
            m_oldTable[i] = EMPTY;
            m_oldNumDeleted++;
        }
    }

    //deletes old table
    if (m_oldNumDeleted == m_oldSize){
        delete [] m_oldTable;
        m_oldTable = nullptr;
    }
}


bool VDetect::remove(Virus virus){
    //checks if the virus is valid
    if (virus.m_id < MINID || virus.m_id > MAXID){
        return false;
    }
    //creates i and sets index to the hash mod capacity
    int i = 0;
    int index = m_hash(virus.m_key) % m_currentCap;

    //iterates while the table at index is not empty and while i < current capacity
    while (!(m_currentTable[index] == EMPTY) && (i <= m_currentCap)){
         if (m_currProbing == NONE){
            return false;
        }

         //calculates the index for quadratic
        else if (m_currProbing == QUADRATIC){
            index = ((m_hash(virus.m_key) % m_currentCap) + (i * i)) % m_currentCap;
        }

        //calculates the index for doubleHash
        else if (m_currProbing == DOUBLEHASH){
            index = ((m_hash(virus.m_key) % m_currentCap) + i * (11 - (m_hash(virus.m_key) % 11))) % m_currentCap;
        }

        //checks if the virus at the index equals the virus to remove and removes it then decrements size
        if (m_currentTable[index] == virus){
            m_currentTable[index] = DELETED;
            m_currentSize--;
            return true;
        }
        i++;
    }

    //checks if the deleted ratio is greater than 0.8 or old table is null
    float deletedCapacity = deletedRatio();
    if ((deletedCapacity > 0.8) or (m_oldTable != nullptr)){
        rehash();
    }

    return false;
}


Virus VDetect::getVirus(string key, int id) const{
    //creates i and sets index to the hash mod capacity
    int i = 0;
    int index = m_hash(key) % m_currentCap;
    Virus emptyVirus("", 0);

    //iterates while the table at index is not empty and while i < current capacity
    while (!(m_currentTable[index] == EMPTY) && (i <= m_currentCap)){
        if (m_currProbing == NONE){
            return emptyVirus;
        }

        //calculates the index for quadratic
        else if (m_currProbing == QUADRATIC){
            index = ((m_hash(key) % m_currentCap) + (i * i)) % m_currentCap;
        }

        //calculates the index for doubleHash
        else if (m_currProbing == DOUBLEHASH){
            index = ((m_hash(key) % m_currentCap) + i * (11 - (m_hash(key) % 11))) % m_currentCap;
        }

        //checks if the index and key matches and then returns the virus at the index
        if (m_currentTable[index].m_id == id && m_currentTable[index].m_key == key ){
            return m_currentTable[index];
        }
        i++;
    }

    if (m_oldTable != nullptr){
        int j = 0;
        int oldIndex = m_hash(key) % m_oldCap;

        //iterates while the table at index is not empty and while j < current capacity
        while (!(m_oldTable[oldIndex] == EMPTY) && (j <= m_oldCap)) {
            if (m_oldProbing == NONE) {
                return emptyVirus;
            }

            //calculates the index for quadratic
            else if (m_oldProbing == QUADRATIC) {
                oldIndex = (m_hash(key) % m_oldCap) + (j * j) % m_oldCap;
            }

            //calculates the index for doubleHash
            else if (m_oldProbing == DOUBLEHASH) {
                oldIndex = ((m_hash(key) % m_oldCap) + j * (11 - (m_hash(key) % 11))) % m_oldCap;
            }

            //checks if the index and key matches and then returns the virus at the index
            if (m_oldTable[oldIndex].m_id == id && m_oldTable[oldIndex].m_key == key) {
                return m_oldTable[oldIndex];
            }
            j++;
        }
    }
    return emptyVirus;
}


float VDetect::lambda() const {
    //calculates the load factor
    float loadFactor = (float)m_currentSize / (float)m_currentCap;
    return loadFactor;
}

float VDetect::deletedRatio() const {
    //calculates the deleted ratio
    float new_ans = m_currNumDeleted / m_currentCap;
    return new_ans;
}

void VDetect::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool VDetect::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int VDetect::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) {
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0)
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

ostream& operator<<(ostream& sout, const Virus &virus ) {
    if (!virus.m_key.empty())
        sout << virus.m_key << " (ID " << virus.m_id << ")";
    else
        sout << "";
    return sout;
}

bool operator==(const Virus& lhs, const Virus& rhs){
    return ((lhs.m_key == rhs.m_key) && (lhs.m_id == rhs.m_id));
}
