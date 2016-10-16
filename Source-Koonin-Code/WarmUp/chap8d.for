       SUBROUTINE LOGDST(X)
       INTEGER SEED
       DATA SEED /98347927/

       X=-LOG(1.-RAN(SEED))
       RETURN
       END
