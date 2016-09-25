#include "Event.h"

void EventHandler::Init(Field* aField){
    mField = aField;
    mMolecule = 0;
    mAtom = 0;

    nIter = 0;
    time = 0;
    timedelta = 0;
}
