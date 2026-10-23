#include "SimTKcommon.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

class ValueFish;
class PointerFish;

static const std::string CHK_STR{"I am a fish"};

class FishImpl {
public:
    FishImpl() :
        m_s("I am a fish")
    {}

    const std::string & getS() const { return m_s; }

private:
    std::string m_s;
};

class ValueFishRep : public SimTK::PIMPLImplementation<ValueFish, ValueFishRep>, public FishImpl {
public:
    ValueFishRep() = default;
    ValueFishRep* clone() const { return new ValueFishRep(); }
    int refCount() const { return getHandleCount(); }
};

class ValueFish : public SimTK::PIMPLHandle<ValueFish, ValueFishRep, false> {
public:
    ValueFish() :
        HandleBase(new ValueFishRep())
    {}

    const std::string & getS() const { return getImpl().getS(); }
    bool isMovedAway() const { return isEmptyHandle(); }
    int refCount() const { return getImpl().refCount(); }
};

class PointerFishRep : public SimTK::PIMPLImplementation<PointerFish, PointerFishRep>, public FishImpl {
public:
    PointerFishRep() = default;
    PointerFishRep* clone() const { return new PointerFishRep(); }
    int refCount() const { return getHandleCount(); }
};

class PointerFish : public SimTK::PIMPLHandle<PointerFish, PointerFishRep, true> {
public:
    PointerFish() :
        HandleBase(new PointerFishRep())
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

    SimTK_TEST(fishCopy.isMovedAway() == false);
    SimTK_TEST(fishCCopy.isMovedAway() == false);
    SimTK_TEST(fishMove.isMovedAway() == false);
    SimTK_TEST(fish.isMovedAway() == true);

    SimTK_TEST(fishCopy.refCount() == 1);
    SimTK_TEST(fishCopy.getS() == CHK_STR);

    SimTK_TEST(fishCCopy.refCount() == 1);
    SimTK_TEST(fishCCopy.getS() == CHK_STR);

    SimTK_TEST(fishMove.refCount() == 1);
    SimTK_TEST(fishMove.getS() == CHK_STR);

    ValueFish fishCMove{std::move(fishCopy)};
    SimTK_TEST(fishCopy.isMovedAway() == true);
}

void testPointerPIMPL()
{
    PointerFish pFish;
    PointerFish pFishCopy = pFish;
    PointerFish pFishCCopy{pFish};

    SimTK_TEST(pFish.isMovedAway() == false);
    SimTK_TEST(pFishCopy.isMovedAway() == false);
    SimTK_TEST(pFishCCopy.isMovedAway() == false);

    SimTK_TEST(pFish.refCount() == 3);
    SimTK_TEST(pFish.getS() == CHK_STR);

    SimTK_TEST(pFishCopy.refCount() == 3);
    SimTK_TEST(pFishCopy.getS() == CHK_STR);

    SimTK_TEST(pFishCCopy.refCount() == 3);
    SimTK_TEST(pFishCCopy.getS() == CHK_STR);

    PointerFish pFishMove = std::move(pFish);
    SimTK_TEST(pFish.isMovedAway() == true);

    SimTK_TEST(pFishCopy.refCount() == 3);
    SimTK_TEST(pFishCopy.getS() == CHK_STR);

    SimTK_TEST(pFishCCopy.refCount() == 3);
    SimTK_TEST(pFishCCopy.getS() == CHK_STR);

    SimTK_TEST(pFishMove.refCount() == 3);
    SimTK_TEST(pFishMove.getS() == CHK_STR);

    PointerFish pFishCMove(std::move(pFishCopy));
    SimTK_TEST(pFishCopy.isMovedAway() == true);

    SimTK_TEST(pFishCCopy.refCount() == 3);
    SimTK_TEST(pFishCCopy.getS() == CHK_STR);

    SimTK_TEST(pFishMove.refCount() == 3);
    SimTK_TEST(pFishMove.getS() == CHK_STR);

    SimTK_TEST(pFishCMove.refCount() == 3);
    SimTK_TEST(pFishCMove.getS() == CHK_STR);
}

int main()
{
    testValuePIMPL();
    testPointerPIMPL();

    return 0;
}
