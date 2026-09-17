// vim: expandtab ts=4 sts=4 sw=4 :

#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

class ValueFish;
class PointerFish;

static const std::string CHK_STR{"I am a fish"};

#define CHKVAL(fish) ASSERT(fish.getS() == CHK_STR)

class FishImpl {
public:
    FishImpl() :
        m_s{"I am a fish"}
    {}

    const std::string & getS() const { return m_s; }

private:
    std::string m_s;
};

class ValueFishRep : public SimTK::PIMPLImplementation<ValueFish, ValueFishRep>, public FishImpl {
public:
    ValueFishRep* clone() const { return new ValueFishRep(); }
    int refCount() const { return getHandleCount(); }
};

class ValueFish : public SimTK::PIMPLHandle<ValueFish, ValueFishRep, false> {
public:
    ValueFish() :
        HandleBase{ new ValueFishRep{} }
    {}

    const std::string & getS() const { return getImpl().getS(); }
    bool isMovedAway() const { return isEmptyHandle(); }
    int refCount() const { return getImpl().refCount(); }
};

class PointerFishRep : public SimTK::PIMPLImplementation<PointerFish, PointerFishRep>, public FishImpl {
public:
    PointerFishRep* clone() const { return new PointerFishRep(); }
    int refCount() const { return getHandleCount(); }
};

class PointerFish : public SimTK::PIMPLHandle<PointerFish, PointerFishRep, true> {
public:
    PointerFish() :
        HandleBase{ new PointerFishRep{} }
    {}

    const std::string & getS() const { return getImpl().getS(); }
    bool isMovedAway() const { return isEmptyHandle(); }
    int refCount() const { return getImpl().refCount(); }
};

void testValuePIMPL()
{
    ValueFish fish;
    ValueFish fishCopy = fish;
    ValueFish fishCCopy{fish};
    ValueFish fishMove = std::move(fish);

    ASSERT(fishCopy.isMovedAway() == false)
    ASSERT(fishCCopy.isMovedAway() == false)
    ASSERT(fishMove.isMovedAway() == false)
    ASSERT(fish.isMovedAway() == true)
    ASSERT(fishCopy.refCount() == 1); CHKVAL(fishCopy)
    ASSERT(fishCCopy.refCount() == 1); CHKVAL(fishCCopy)
    ASSERT(fishMove.refCount() == 1); CHKVAL(fishMove)

    ValueFish fishCMove{std::move(fishCopy)};
    ASSERT(fishCopy.isMovedAway() == true)
}

void testPointerPIMPL()
{
    PointerFish pFish;
    PointerFish pFishCopy = pFish;
    PointerFish pFishCCopy{pFish};

    ASSERT(pFish.isMovedAway() == false)
    ASSERT(pFishCopy.isMovedAway() == false)
    ASSERT(pFishCCopy.isMovedAway() == false)
    ASSERT(pFish.refCount() == 3); CHKVAL(pFish)
    ASSERT(pFishCopy.refCount() == 3); CHKVAL(pFishCopy)
    ASSERT(pFishCCopy.refCount() == 3); CHKVAL(pFishCCopy)

    PointerFish pFishMove = std::move(pFish);
    ASSERT(pFish.isMovedAway() == true);
    ASSERT(pFishCopy.refCount() == 3); CHKVAL(pFishCopy)
    ASSERT(pFishCCopy.refCount() == 3); CHKVAL(pFishCCopy);
    ASSERT(pFishMove.refCount() == 3); CHKVAL(pFishMove);

    PointerFish pFishCMove(std::move(pFishCopy));
    ASSERT(pFishCopy.isMovedAway() == true)
    ASSERT(pFishCCopy.refCount() == 3); CHKVAL(pFishCCopy)
    ASSERT(pFishMove.refCount() == 3); CHKVAL(pFishMove);
    ASSERT(pFishCMove.refCount() == 3); CHKVAL(pFishCMove);
}

int main()
{
    testValuePIMPL();
    testPointerPIMPL();

    return 0;
}
