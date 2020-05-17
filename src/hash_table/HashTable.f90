MODULE HashTable

  USE big_integer_module

CONTAINS

  SUBROUTINE key2key(key1, key2, reaction, reactionkey, rkeysign, pd)
    USE big_integer_module
    IMPLICIT NONE

    INTEGER reaction, pd, rkeysign(pd)
    TYPE (big_integer) :: key1, key2, reactionkey(pd)

    IF (rkeysign(reaction).GT.0) THEN
       key2 = key1 + reactionkey(reaction)
    ELSE
       key2 = key1 - reactionkey(reaction)
    ENDIF
  END SUBROUTINE  key2key

  SUBROUTINE key2keybw(key1, key2, reaction, reactionkey, rkeysign, pd)
    USE big_integer_module
    IMPLICIT NONE

    INTEGER reaction, pd, rkeysign(pd)
    TYPE (big_integer) :: key1, key2, reactionkey(pd)

    IF (rkeysign(reaction)>0) THEN
       IF (key1<reactionkey(reaction)) THEN
          key2 = 0
       ELSE
          key2 = key1 - reactionkey(reaction)
       ENDIF
    ELSE
       key2 = key1 + reactionkey(reaction)
    ENDIF
  END SUBROUTINE key2keybw

  SUBROUTINE state2key(KEY, state, sd, B)
    !---  returns the enumeration of the state
    USE big_integer_module
    IMPLICIT NONE
    INTEGER sd, state(sd), B, k
    TYPE (big_integer) :: j, KEY, C

    C = big(B)
    j = 2
    DO k = 1, sd
       j = j + state(k) * (C + 1)**(k - 1)
    ENDDO
    DO k = 1, sd
       !        check state is
       !        1.   nonnegative;  and
       !        2.   inside maximum allowed;
       !        (otherwise set j=0 as flag)
       IF (state(k).LT.0 .OR. state(k).GT.B) j = 0
    ENDDO
    KEY = j
  END SUBROUTINE state2key
  !----------------------------------------------------------------------|
  SUBROUTINE HASH(KEY, MODE, KA, FOUND, &
       KTLEN, KEYTAB, KVTAB)
    !     Adapted for our cme-solver.
    !     Uses big_integer_module
    !     Should take advantage of the fact that DELETE is not used.
    !     Changes added:
    !     - A constant DELKEY for deleted key (since big_integer does not
    !     support negative numbers).
    !     - Implicit none
    !===============BRENT'S HASHING ALGORITHM.
    !               SEE CACM 16(2):105-109.
    !
    !       THIS ROUTINE WILL LOOK UP AND INSERT OR DELETE
    !       INTEGER KEYS IN THE TABLE KEYTAB.  VALUES ASSOCIATED
    !       WITH THE KEYS MAY BE STORED IN THE TABLE KVTAB. THE
    !       ROUTINE IS DESIGNED TO BE EFFICIENT, EVEN IF THE TABLE
    !       IS NEARLY FULL, PROVIDED MOST ENTRIES ARE LOOKED UP
    !       SEVERAL TIMES.
    !  PARAMETERS:
    !  /BRENT/ COMMON BLOCK DESCRIBES:
    !          LEN:    THE TABLE LENGTH.  MUST BE GIVEN AS A PRIME NUMBER
    !                  BEFORE HASH IS CALLED. NO CHECK IS MADE FOR PRIMALITY
    !
    !          KEYTAB: IS GIVEN AS A LEN-VECTOR FOR STORING THE KEYS. IT
    !                  MUST BE INITIALIZED TO ALL ZEROS BEFORE THE FIRST
    !                  CALL TO HASH.
    !
    !          KVTAB:  IS USED BY BOTH HASH AND THE CALLING PROGRAM AS
    !                  A LEN-VECOTR OF VALUES ASSOCIATED WITH KEYS IN
    !                  THE HASH TABLE, KEYTAB.  KVTAB COULD ALSO BE A
    !                  POINTER VECTOR INITIALIZED TO THE INTEGERS 1,2,...LEN
    !                  IF MULTIPLE VALUES ARE ASSOCIATED WITH EACH KEY.
    !  (KEY,MODE,KA,FOUND):
    !          KEY:    IS GIVEN AS A POSITIVE INTEGER KEY TO BE ENTERED,
    !                  DELETED, OR PROBED FOR.
    !                  THE VALUES 0 AND -1 ARE RESERVED FOR EMPTY SPACE
    !                  OR DELETED ITEM RESPECTIVELY.
    !
    !          MODE:   IS GIVEN AS AN INDICATOR AS FOLLOWS...
    !          MODE=1: MEANS LOOK UP ONLY.  IF THE KEY IS IN KEYTAB, THEN
    !                  KA IS RETURNED AS THE SUBSCRIPT POINTER; OTHERWISE
    !                  KA IS RETURNED 0.  THE VALUE ASSOCIATED WITH
    !                  KEYTAB(KA) IS KVTAB(KA) WHEN KA IS NOT ZERO.
    !          MODE=2: MEANS LOOK UP AND ENTER.  THIS IS USED TO BUILD
    !                  THE TABLE.  OTHERWISE THE KEY IS ENTERED AT
    !                  KEYTAB(KA).  KA IS RETURNED AS THE POINTER SUBSCRIPT
    !                  IN EITHER CASE, OR KA=0 IF AN ENTRY IS ATTEMPTED AND
    !                  THE TABLE IS FULL OR KA=0 IF KEY=0.  UPON RETURN,
    !                  WHEN FOUND=.FALSE., THE CALLING PROGRAM MUST ENTER
    !                  THE VALUE ASSOCIATED WITH KEY WHEN MODE=2.
    !          MODE=3: MEANS LOOK UP AND DELETE.  IF THE KEY IS IN THE TABLE
    !                  IT IS DELETED AND ITS FORMER ADDRESS KA IS RETURNED.
    !                  IF THE KEY IS NOT THERE, KA=0 IS RETURNED.  (THE
    !                  SUBROUTINE CAN BE SIMPLIFIED CONSIDERABLY IF KEYS
    !                  NEVER DELETED).
    !
    !          KA:     IS RETUNED AS THE SUBSCRIPT OF KEYTAB SUCH THAT
    !                  KEY=KEYTAB(KA) OR KA IS RETURNED AS ZERO IF KEY=0
    !                  OR KEY IS NOT IN KEYTAB.
    !
    !          FOUND:  IS A LOGICAL FLAG RETURNED .TRUE. IF THE KEY WAS
    !                  FOUND IN KEYTAB; ELSE FOUND IS RETURNED .FALSE. .
    !
    !          KTLEN: LENGTH OF KEY TABLE

    USE big_integer_module

    IMPLICIT NONE

    INTEGER IQ, IR, IC, LEN2, IA, IX, JQ, JR
    LOGICAL FOUND, DEL
    !-----------THE COMMON STATEMENT  BELOW MUST BE CHANGED TO MEET
    !           THE PROBLEM SIZE.
    INTEGER KTLEN, MODE, KA
    INTEGER KVTAB(KTLEN)
    TYPE (big_integer) :: KEY, KEYTAB(KTLEN), IS, KT, DELKEY
    !      INTEGER KEY, KEYTAB(KTLEN)

    DELKEY = 1

    LEN2 = KTLEN - 2
    IC = -1
    !--------COMPUTE ADDRESS OF FIRST PROBE(IR) AND INCREMENT(IQ).
    !        ANY INDEPENDENT PSEUDO-RANDOM FUNCTIONS OF THE KEY MAY BE USED,
    !        PROVIDED 0<IQ<LEN AND 0<IR<=LEN .
    IQ = MODULO(KEY, LEN2) + 1
    IR = MODULO(KEY, KTLEN) + 1
    !      IQ=MOD(IABS(KEY),LEN2)+1
    !      IR=MOD(IABS(KEY),KTLEN)+1
    KA = IR
    !-----------LOOK IN THE TABLE.
20  KT = KEYTAB(KA)
    !----------CHECK FOR AN EMPTY SPACE, A DELETE ENTRY, OR A MATCH.
    IF(KT.EQ.0) GO TO 30
    IF(KT.EQ.DELKEY) GO TO 40
    IF(KT.EQ.KEY) GO TO 60
    IC = IC + 1
    !---------COMPUTE ADDRESS OF NEXT PROBE.
    KA = KA + IQ
    IF(KA.GT.KTLEN)KA = KA - KTLEN
    !----------SEE IF WHOLE TABLE HAS BEEN SEARCHED.
    IF(KA.NE.IR) GO TO 20
    !---------THE KEY IS NOT IN THE TABLE.
30  FOUND = .FALSE.
    !---------RETURN WITH KA=0 UNLES AN ENTRY HAS TO BE MADE.
    IF((MODE.EQ.2).AND.(IC.LE.LEN2) .AND.(KEY.NE.0).AND.&
         (KEY.NE.DELKEY))&
         GO TO 70
    KA = 0
    RETURN
    !-----------A DELETED ENTRY HAS BEEN FOUND.
40  IA = KA
    !----------COMPUTE ADDRESS OF NEXT PROBE.
50  IA = IA + IQ
    IF(IA.GT.KTLEN) IA = IA - KTLEN
    IS = KEYTAB(IA)
    !----------CHECK FOR AN EMPTY SPACE OR A COMPLETE SCAN OF THE TABLE.
    IF((IS.EQ.0) .OR.(IA.EQ.IR)) GO TO 30
    !------------CHECK FOR A MISMATCH OR A DELETED ENTRY.
    IF((IS.NE.KEY).OR.(IS.EQ.DELKEY)) GO TO 50
    !----------KEY FOUND.  MOVE IT AND THE ASSOCIAATED VALUE TO SAVE PROBES
    !          ON THE NEXT SEARCH FOR THE SAME KEY.
    KVTAB(KA) = KVTAB(IA)
    KEYTAB(KA) = IS
    KEYTAB(IA) = DELKEY
    !----------THE KEY IS IN THE TABLE.
60  FOUND = .TRUE.
    !----------DELETE IT IF MODE=3
    IF(MODE.EQ.3) KEYTAB(KA) = DELKEY
    RETURN
    !-----------LOOK FOR THE BEST WAY TO MAKE AN ENTRY.
70  IF(IC.LE.0) GO TO 120
    !----------SET DEL IF A DELETED ENTRY HAS BEEN FOUND.
    DEL = KT.NE.0
    IA = KA
    IS = 0
    !-----------COMPUTE THE MAXIMUM LENGTH TO SEARCH ALONG CURRENT CHAIN.
80  IX = IC - IS
    !-----------COMPUTE INCREMENT JQ FOR CURRENT CHAIN.
    !      JQ=MOD(IABS(KEYTAB(IR)),LEN2)+1
    JQ = MODULO(KEYTAB(IR), LEN2) + 1
    JR = IR
    !----------LOOK ALONG THE CHAIN.
90  JR = JR + JQ
    IF(JR.GT.KTLEN) JR = JR - KTLEN
    KT = KEYTAB(JR)
    !----------CHECK FOR A HOLE (AN EMPTY SPACE OR DELETED ENTRY).
    IF((KT.EQ.0).OR.(KT.EQ.DELKEY)) GO TO 100
    IX = IX - 1
    IF(IX.GT.0) GO TO 90
    GO TO 110
    !----------SKIP IF THIS IS AN EMPTY SPACE AND A DELETED ENTRY HAS
    !          ALREADY BEEN FOUND.
100 IF(DEL.AND.(KT.EQ.0)) GO TO 110
    !------------CHECK FOR A DELETED ENTRY.
    IF(KT.NE.0) DEL = .TRUE.
    !-----------SAVE LOCATION OF HOLE.
    IA = JR
    KA = IR
    IC = IC - IX
    !------------MOVE DOWN TO THE NEXT CHAIN.
110 IS = IS + 1
    IR = IR + IQ
    IF(IR.GT.KTLEN) IR = IR - KTLEN
    !---------GO BACK IF A BETTER HOLE MIGHT STILL BE FOUND.
    IF(IC.GT.IS) GO TO 80
    !---------SKIP IF THERE IS NOTHING TO MOVE.
    IF(IA.EQ.KA) GO TO 120
    !---------MOVE AN OLD ENTRY AND ITS ASSOCIATED VALUE TO MAKE ROOM FOR
    !         THE NEW ENTRY.
    KVTAB(IA) = KVTAB(KA)
    KEYTAB(IA) = KEYTAB(KA)
    !---------ENTER THE NEW KEY, BUT NOT ITS ASSICIATED VALUE.
120 KEYTAB(KA) = KEY
    RETURN
  END SUBROUTINE HASH

END MODULE HashTable