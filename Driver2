#include "vdetect.h"
#include <random>
#include <vector>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};
class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

class Tester{
public:
    bool testInsertionOperation();
    bool testFindErrorCase();
    bool testFindNonCollidingKeys();
    bool testFindCollidingKeys();
    bool testRemoveNonColliding();
    bool testRemoveColliding();
    bool testRehashInsertion();
    bool testRehashLoadFactor();
    bool testRehashRemoval();
    bool testRehashDeleteRatio();
private:
};

int main(){
    Tester Tester1;

    if (Tester1.testInsertionOperation()){
        cout << "Test Insert Normal returned true, works correctly for normal case" << endl;
    }
    else{
        cout << "Test Insert normal returned false" << endl;
    }

    if (Tester1.testFindErrorCase()){
        cout << "Test Find Error Case returned true, works correctly" << endl;
    }
    else{
        cout << "Test Find Error Case returned false" << endl;
    }

    if (Tester1.testFindNonCollidingKeys()){
        cout << "Test Find Error Non Colliding keys returned true, works correctly" << endl;
    }
    else{
        cout << "Test Find Non colliding keys returned false" << endl;
    }

    if (Tester1.testFindCollidingKeys()){
        cout << "Test Find Error Colliding keys returned true, works correctly" << endl;
    }
    else{
        cout << "Test Find colliding keys returned false" << endl;
    }

    if (Tester1.testRemoveNonColliding()){
        cout << "Test Remove Non Colliding keys returned true, works correctly" << endl;
    }
    else{
        cout << "Test Remove Non colliding keys returned false" << endl;
    }

    if (Tester1.testRemoveColliding()){
        cout << "Test Remove Colliding keys returned true, works correctly" << endl;
    }
    else{
        cout << "Test Remove colliding keys returned false" << endl;
    }

    if (Tester1.testRehashInsertion()){
        cout << "Test Rehash Insertion returned true, works correctly" << endl;
    }
    else{
        cout << "Test Rehash Insertion returned false" << endl;
    }

    if (Tester1.testRehashLoadFactor()){
        cout << "Test Rehash Load Factor returned true, works correctly" << endl;
    }
    else{
        cout << "Test Rehash Load Factor returned false" << endl;
    }

    if (Tester1.testRehashRemoval()){
        cout << "Test Rehash Removal returned true, works correctly" << endl;
    }
    else{
        cout << "Test Rehash Removal returned false" << endl;
    }

    if (Tester1.testRehashDeleteRatio()){
        cout << "Test Rehash Delete Ratio returned true, works correctly" << endl;
    }
    else{
        cout << "Test Rehash Delete Ratio returned false" << endl;
    }

}


bool Tester::testInsertionOperation(){
    bool ans = true;
    bool result = false;

    vector<Virus> dataList;
    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 54;

    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }
    // dumping the data points to the standard output
    //vdetect.dump();

    // checking whether all data are inserted
    for (vector<Virus>::iterator it = dataList.begin(); it != dataList.end(); it++){
        ans = ans && (*it == vdetect.getVirus((*it).getKey(), (*it).getID()));
    }

    //check whether the size is updated correctly;
    ans = ans && (vdetect.m_currentSize == size);

    if (ans){
        result = true;
        return result;
    }

    return result;
}


bool Tester::testFindErrorCase(){
    bool ans = false;
    Virus invalidVirus("AEIOU", 1543);

    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 54;

    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }

    //finding a virus that does not exist in database
    Virus tempVirus = vdetect.getVirus(invalidVirus.getKey(), invalidVirus.getID());
    //checks if the function returns empty string and -1
    if ( tempVirus.getKey() == "" && tempVirus.getID() == 0){
        ans = true;
    }

    return ans;
}


bool Tester::testFindNonCollidingKeys(){
    bool ans = false;
    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 10;

    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }
    //creates temp object and inserts into the virus
    Virus tempObj = Virus(sequencer(5, 1), 3456);
    vdetect.insert(tempObj);
    //stores the key and id
    string tempKey = tempObj.getKey();
    int tempId = tempObj.getID();
    //finds the virus and checks if it equals the tempObj
    Virus currVirus = vdetect.getVirus(tempKey, tempId);
    if (tempObj == currVirus){
        ans = true;
    }
    return ans;
}


bool Tester::testFindCollidingKeys(){
    bool ans = false;
    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 5;

    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }
    //Creates temp viruses
    Virus tempVirus1 = Virus("ABCDE", 3456);
    Virus tempVirus2 = Virus("JOLO", 4567);
    Virus tempVirus3 = Virus("EXRP", 8910);
    Virus tempVirus4 = Virus("ANMN", 1986);
    Virus tempVirus5 = Virus("AQR2", 4848);
    //inserts into the vdetect
    vdetect.insert(tempVirus1);
    vdetect.insert(tempVirus2);
    vdetect.insert(tempVirus3);
    vdetect.insert(tempVirus4);
    vdetect.insert(tempVirus5);

    //creates temp viruses with different keys but same ID
    Virus tempVirus11 = Virus("AfgE", 3456);
    Virus tempVirus22 = Virus("gkjR", 4567);
    Virus tempVirus33 = Virus("klRc", 8910);
    Virus tempVirus44 = Virus("HUmn", 4098);
    Virus tempVirus55 = Virus("dfHL", 4848);
    //inserts into the vdetect
    vdetect.insert(tempVirus11);
    vdetect.insert(tempVirus22);
    vdetect.insert(tempVirus33);
    vdetect.insert(tempVirus44);
    vdetect.insert(tempVirus55);

    string tempKey = tempVirus11.getKey();
    int tempId = tempVirus11.getID();

    //checks if the colliding virus was inserted
    Virus currVirus = vdetect.getVirus(tempKey, tempId);
    if (tempVirus11 == currVirus){
        ans = true;
    }

    return ans;
}


bool Tester::testRemoveNonColliding(){
    bool ans = false;

    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 10;

    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }
    //creates temp object and inserts into the virus
    Virus tempObj = Virus(sequencer(5, 1), 3456);
    vdetect.insert(tempObj);
    //stores the key and id
    string tempKey = tempObj.getKey();
    int tempId = tempObj.getID();
    //remove the virus
    vdetect.remove(tempObj);

    //finds the virus and checks if it is not found then it is deleted
    Virus currVirus = vdetect.getVirus(tempKey, tempId);
    if ((currVirus.getKey() == "") && (currVirus.getID() == 0)){
        ans = true;
    }

    return ans;
}


bool Tester::testRemoveColliding(){

    bool ans = false;
    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 5;

    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }
    //Creates temp viruses
    Virus tempVirus1 = Virus("ABCDE", 3456);
    Virus tempVirus2 = Virus("JOLO", 4567);
    Virus tempVirus3 = Virus("EXRP", 8910);
    Virus tempVirus4 = Virus("ANMN", 1986);
    Virus tempVirus5 = Virus("AQR2", 4848);
    //inserts into the vdetect
    vdetect.insert(tempVirus1);
    vdetect.insert(tempVirus2);
    vdetect.insert(tempVirus3);
    vdetect.insert(tempVirus4);
    vdetect.insert(tempVirus5);

    //creates temp viruses with different keys but same ID
    Virus tempVirus11 = Virus("AfgE", 3456);
    Virus tempVirus22 = Virus("gkjR", 4567);
    Virus tempVirus33 = Virus("klRc", 8910);
    Virus tempVirus44 = Virus("HUmn", 4098);
    Virus tempVirus55 = Virus("dfHL", 4848);
    //inserts into the vdetect
    vdetect.insert(tempVirus11);
    vdetect.insert(tempVirus22);
    vdetect.insert(tempVirus33);
    vdetect.insert(tempVirus44);
    vdetect.insert(tempVirus55);

    string tempKey = tempVirus11.getKey();
    int tempId = tempVirus11.getID();
    //remove the virus
    vdetect.remove(tempVirus11);
    //finds the virus and checks if it is not found then it is deleted
    Virus currVirus = vdetect.getVirus(tempKey, tempId);
    if ((currVirus.getKey() == "") && (currVirus.getID() == 0)){
        ans = true;
    }
    return ans;
}


bool Tester::testRehashInsertion(){
    bool ans = false;

    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 51;
    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }

    //checks if the old table is not null that means the rehash was successful
    if(vdetect.m_oldTable != nullptr){
        ans = true;
    }

    return ans;
}


bool Tester::testRehashLoadFactor(){
    bool ans = false;

    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 54;
    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }

    if(vdetect.lambda() <= 0.5){
        ans = true;
    }
    return ans;
}


bool Tester::testRehashRemoval(){
    bool ans = false;

    vector<Virus> dataList;
    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 50;
    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);

        // vector to store virus
        dataList.push_back(dataObj);
    }

    for (unsigned int i = 0; i < dataList.size() / 3; i++){
        vdetect.remove(dataList[i]);
    }
    //vdetect.dump();
    if(vdetect.m_oldTable == nullptr){
        ans = true;
    }
    return ans;
}


bool Tester::testRehashDeleteRatio(){
    bool ans = false;

    Random RndID(MINID,MAXID);
    VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
    int size = 54;
    for (int i = 0; i < size; i++){
        // generating random data
        Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
        // inserting data in to the VDetect object
        vdetect.insert(dataObj);
    }

    if(vdetect.deletedRatio() <= 0.8){
        ans = true;
    }
    return ans;
}


unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;  // magic number from textbook
    for ( int i = 0 ; i < str.length(); i++)
        val = val * thirtyThree + str[i] ;
    return val ;
}

string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}
