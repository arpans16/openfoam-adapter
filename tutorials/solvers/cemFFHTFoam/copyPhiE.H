        double time2 = runTime.value();
        if (time1 == time2)
        {
            count++;
        }
        else
        {
            count = 1;
            time1 = time2;
        }
        if (count == 1)
        {
            phiEold == phiE.oldTime();
        }
        else
        {
            phiEold = phiE.oldTime();
        }
