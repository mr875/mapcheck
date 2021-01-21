from .types import Types
from datetime import datetime
import sys
import os
import re
from utils.queryfile import QueryFile, NormFile

class ProcFile(Types):

    def __init__(self,fname,tabname,omics,br):
        self.inp = NormFile(fname)
        self.ts = datetime.now().strftime("%d%b_%Hhr")
        self.omics = omics
        self.br = br
        self.tabname = tabname
        brk=10
        print("line count",self.inp.row_count)
        self.dbact_om = open('dbact_om_' + self.ts + '.sql',"w")
        self.dbact_br = open('dbact_br_' + self.ts + '.sql',"w")
        if 'map_new_alt_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.newaltrs(brk=brk)
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_alt_rs.sh',self.inp.bfile,'out_map_new_alt_rs.txt'))
        self.dbact_om.close()
        self.dbact_br.close()

    def grabrsinp(self,col):
        rs = re.search("(?:\()(rs[0-9]+)(?:\))", col)
        if rs:
            return rs.group(1)
        else:
            raise Exception("grabrsinp(col): can't find rs id in parenthesis: %s" % (col))

    def mergecheck(self,col):
        rside = re.search("(?:[mM]erges?:)(.+)",col)
        rside = rside.group(1)
        if 'rs' not in rside:
            return {'withdrawn':False,'merges':[]}
        if 'withdrawn' in rside:
            return {'withdrawn':True}
        mg = re.findall("(rs[0-9]+)",rside)
        if (len(mg) % 2) != 0:
            raise Exception("mergecheck(col): expecting paired (even) rsids in merge list: ", col)
        return {'withdrawn':False,'merges':mg}

    def frominto(self,evenarr,bfor,aftr):
        #print(evenarr)
        nxt = []
        correct = False
        after_match = False
        for ind,rs in enumerate(evenarr):
            ind += 1
            if (ind % 2) != 0:
                if rs == bfor:
                    correct = True
                    nxt.append(evenarr[ind])
        if correct:
            if aftr in nxt:
                after_match = True
        return {'correct':correct,'after_match':after_match,'alt_after':nxt}

    def getb38(self,container):
        rslt = re.search('(?:b38=)(\S+)',container)
        return rslt.group(1) # will need handling: AttributeError: 'NoneType' object has no attribute 'group'

    def add_alt(self,alt,main,ds):
        pass

    def mtab_change_id(self,xsting,chngto):
        pass

    def mtab_get_where_string(self,anid):
        if 'ukbbaffy_v2_1_map' in self.tabname:
            rscol = 'chipid'
        else:
            rscol = 'dbsnpid'
        where = ' WHERE snp = %s'
        q = 'SELECT snp,dbsnpid,chr FROM ' + self.tabname + where
        val = (anid,)
        self.br.execute(q,val)
        res = self.br.fetchall()
        if self.br.rowcount == 0 and 'rs' in anid: # or len(res)
            where = ' WHERE ' + rscol + ' REGEXP %s'
            q = 'SELECT snp,dbsnpid,chr FROM ' + self.tabname + where
            val = (anid+'[[:>:]]',)
            self.br.execute(q,val)
            res = self.br.fetchall()
        if self.br.rowcount > 0:
            return where,val,rscol

    def mtab_change_pos(self,anid,oldpos,newpos):
        where,val,rscol = self.mtab_get_where_string(anid)
        q = 'UPDATE ' + self.tabname + ' SET chr = %s ' + where + ' AND chr = %s'
        val = (newpos,) + val + (oldpos,)
        #print(q % val)
        #self.br.execute(q,val)
        #res = self.br.fetchall()
        #res = ['|'.join(row) for row in res]

    def swapout_main(self,swin,swout,ds):
        q = 'SELECT uid_datasource FROM consensus WHERE id = %s'
        val = (swout,)
        self.omics.execute(q,val)
        res = self.omics.fetchall()
        res = [rs[0] for rs in res if len(rs) > 0]
        if len(res) != 1:
            raise Exception('failed/unexpected result:' + q % val)
        swout_ds = res[0]
        self.dbact_om.write('UPDATE consensus SET id = \'%s\', uid_datasource = \'%s\' WHERE id = \'%s\';\n' % (swin,ds,swout))
        self.dbact_om.write('UPDATE alt_ids SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE flank SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE positions SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE probes SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE snp_present SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('INSERT INTO alt_ids (id, alt_id, datasource) VALUES (\'%s\',\'%s\',\'%s\');\n' % (swin,swout,swout_ds))
        self.dbact_om.write('--\n')

    def pos_flag(self,mid,chrpos,fl=-5,supp=True):
        vals = (fl,mid,'38',chrpos,-1)
        q = 'UPDATE positions SET chosen = %s WHERE id = %s AND build = %s AND CONCAT(chr,":",pos) = %s AND chosen > %s'
        if supp:
            return
        self.omics.execute(q,vals)

    def addpos(self,mid,chrpos,ds,build='38',supp=True):
        chrm = chrpos.split(':')[0]
        pos = int(chrpos.split(':')[1])
        q = 'INSERT IGNORE INTO positions (id, chr, pos, build, datasource) VALUES (%s,%s,%s,%s,%s)'
        vals = (mid,chrm,pos,build,ds)
        if supp:
            return
        self.omics.execute(q,vals)
    
    def checkbr_pos(self,mid,chrpos=None):
        where,val,rscol = self.mtab_get_where_string(mid)
        q = 'SELECT chr FROM ' + self.tabname + where
        self.br.execute(q,val)
        res = self.br.fetchall()
        res = [s[0] for s in res]
        if not chrpos:
            return [None,res]
        if chrpos in res:
            return [True,res]
        return [False,res]

    def checkom_pos(self,mid,chrpos=None):
        q = 'SELECT chr,pos FROM positions WHERE build = %s AND id = %s'
        val = ('38',mid)
        self.omics.execute(q,val)
        res = self.omics.fetchall()
        res = [s[0] + ':'+ str(s[1]) for s in res]
        if not chrpos:
            return [None,res]
        if chrpos in res:
            return [True,res]
        return [False,res]

    def extra_map(self,newid,linkid,chrpos,datasource,chosen,ds_chrpos=None):
        chrm = chrpos.split(':')[0]
        pos = chrpos.split(':')[1]
        if not ds_chrpos:
            q = "INSERT IGNORE INTO extra_map (newid,linkid,chr,GRCh38_pos,datasource,chosen) VALUES ('%s','%s','%s','%s','%s',%s)"
            vals = (newid,linkid,chrm,pos,datasource,str(chosen))
        else:
            q = "INSERT IGNORE INTO extra_map (newid,linkid,chr,GRCh38_pos,datasource,ds_chrpos,chosen) VALUES ('%s','%s','%s','%s','%s','%s',%s)"
            vals = (newid,linkid,chrm,pos,datasource,ds_chrpos,str(chosen))
        self.extra_map_f.write((q+';\n') % vals)

