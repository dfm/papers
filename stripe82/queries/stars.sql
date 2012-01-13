SELECT
    p.objID,p.ra,p.dec,p.g,p.Err_g,p.u,p.r,p.i,p.z
INTO mydb.stars
FROM Stripe82..PhotoPrimary p
WHERE p.ra BETWEEN 0.000000 AND 10.000000
AND p.dec BETWEEN -2.000000 AND 2.000000
AND p.type = 6 AND p.g BETWEEN 15.000000 AND 20.000000
AND (p.flags &
    (dbo.fPhotoFlags('BRIGHT')+dbo.fPhotoFlags('EDGE')+dbo.fPhotoFlags('BLENDED')
    +dbo.fPhotoFlags('SATURATED')+dbo.fPhotoFlags('NOTCHECKED')
    +dbo.fPhotoFlags('NODEBLEND')+dbo.fPhotoFlags('INTERP_CENTER')
    +dbo.fPhotoFlags('DEBLEND_NOPEAK')+dbo.fPhotoFlags('PEAKCENTER'))) = 0

