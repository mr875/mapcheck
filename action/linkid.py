

dsmap = {'170':['coreexome_map'],'171':['coreexome_map'],'232':['coreexome_map'],'69':['humanexome_map'],'233':['humanexome_map'],'70':['infiniumimmunoarray_map'],'234':['infiniumimmunoarray_map'],'72':['msexome_map'],'73':['omniexpress_map','omniexpress_v2_1_map'],'235':['omniexpress_map','omniexpress_v2_1_map'],'71':['ukbbaffy_map'],'231':['ukbbaffy_v2_1_map'],'999':[]} 
# NOTE ds 999 has no map table equivalent

def uidcol(maptab):
    uidcols = {'ukbbaffy_v2_1_map':'chipid'}
    return uidcols.get(maptab,'snp') # even for ukbbaffy_map, although its snp == v2_1.snp, no chipid

def rscol(maptab):
    rscols = {}
    return rscols.get(maptab,'dbsnpid')

class LinkID:
    
    @classmethod
    def get_flank_ds(cls,mid,omcurs):
        qry = 'SELECT datasource, chosen FROM flank WHERE id = %s'
        omcurs.execute(qry,(mid,))
        ls = omcurs.fetchall()
        return ls

    @classmethod
    def anotds_to_mapt(cls,ds):
        return dsmap.get(ds,ds) # default = ds, because if it's not in the dict then the ds is probably a map table itself

    @classmethod
    def where_col(cls,mid,maptable,brcurs):
        mt_col = uidcol(maptable)
        rs_col = rscol(maptable)
        where = ' WHERE ' + rs_col + ' REGEXP %s'
        val = (mid+'[[:>:]]',)
        colwithid = rs_col
        res = None
        if 'rs' in mid:
            q = 'SELECT ' + mt_col + ',' + rs_col + ',chr FROM ' + maptable + where
            brcurs.execute(q,val)
            res = brcurs.fetchall()
        if not res:
            where = ' WHERE ' + mt_col + ' = %s'
            val = (mid,)
            q = 'SELECT ' + mt_col + ',' + rs_col + ',chr FROM ' + maptable + where
            brcurs.execute(q,val)
            res = brcurs.fetchall()
            colwithid = mt_col
        return [where,val,colwithid,res] # check if val in not None to detect success

    @classmethod
    def get_maptab_where(cls,mid,maptable,brcurs,omcurs=None):
        resarr = cls.where_col(mid,maptable,brcurs)
        if not resarr[3] and omcurs:
            q = 'SELECT alt_id FROM alt_ids WHERE id = %s'
            val = (mid,)
            omcurs.execute(q,val)
            alts = omcurs.fetchall()
            alts = [r[0] for r in alts]
            for alt in alts:
                resarr = cls.where_col(alt,maptable,brcurs)
                if resarr[3]:
                    break
        return resarr
