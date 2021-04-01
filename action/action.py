from .omrs_13 import OmRs_13
from .newrs_07 import NewRs_07
from .newaltrs_05 import NewAltRs_05
from .newposmism_06 import NewPosMisM_06
from .newrsbyalt_08 import NewRsbyAlt_08
from .no38pos_09 import No38Pos_09
from datetime import datetime
import sys
import os
import re
from utils.queryfile import QueryFile, NormFile

class ProcFile(NewAltRs_05,NewPosMisM_06,NewRs_07,NewRsbyAlt_08,No38Pos_09,OmRs_13):

    def __init__(self,fname,tabname,omics,br,reportmode=True):
        self.inp = NormFile(fname)
        self.ts = datetime.now().strftime("%d%b_%Hhr")
        self.omics = omics
        self.br = br
        self.tabname = tabname
        self.reportmode = reportmode
        self.make_extra_map_table()
        brk=0 # set to 0 for whole file # if 20 then line number 19 gets done, line 20 does not
        start=0 # set to 0 for no action. if 10 then line number 10 gets done
        self.dbact_om = open('dbact_om_' + self.ts + '_' + tabname + '.sql',"w")
        self.dbact_br = open('dbact_br_' + self.ts + '_' + tabname + '.sql',"w")
        self.actbr_poschng = set()
        self.actbr_idchng = set()
        self.ignorepos = ['flank_error','pending']
        if 'map_new_alt_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.newaltrs(brk=brk,start=start)
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_alt_rs.sh',self.inp.bfile,'out_map_new_alt_rs.txt'))
        if 'new_pos_mismatch_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.newposmims(brk=brk,start=start)
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_pos_mismatch.txt',self.inp.bfile,'out_new_pos_mismatch_rs.txt'))
        if '_new_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                if 'byalt' not in os.path.basename(self.inp.bfile):
                    self.newrs(brk=brk,start=start)
                else:
                    self.newrsbyalt(brk=brk,start=start)
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_rs.sh',self.inp.bfile,'out_new_rs_rs.txt'))
        if 'no38pos' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.no38pos(brk=brk,start=start)
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_no38pos.txt',self.inp.bfile,'out_no38pos_rs.txt'))
        if 'rsomics' in self.inp.bfile:
            self.omrs(brk=brk,start=start)
        self.dbact_om.close()
        self.dbact_br.close()

    def make_extra_map_table(self):
        if self.reportmode:
            return
        q = 'SELECT COUNT(*) FROM information_schema.tables WHERE table_name = %s'
        val = ('extra_map',)
        self.omics.execute(q,val)
        if self.omics.fetchone()[0] == 0:
            maketable = 'CREATE TABLE IF NOT EXISTS extra_map (newid char(20) NOT NULL,linkid char(38) NOT NULL,chr char(4) DEFAULT NULL,GRCh38_pos int unsigned DEFAULT NULL,flank_seq varchar(1040) DEFAULT NULL,datasource varchar(40) NOT NULL,ds_chrpos varchar(20) DEFAULT NULL,chosen TINYINT DEFAULT 0,PRIMARY KEY (newid,linkid,datasource))'
            self.omics.execute(maketable)

    def stillmain(self,mid):
        q = 'SELECT EXISTS (SELECT 1 FROM consensus WHERE id = %s)'
        val = (mid,)
        self.omics.execute(q,val)
        exists = self.omics.fetchone()[0]
        return exists

    def grabrsinp(self,col):
        rs = re.search("(?:\()(rs[0-9]+)(?:\))", col)
        if rs:
            return rs.group(1)
        else:
            raise Exception("grabrsinp(col): can't find rs id in parenthesis: %s" % (col))

    def mergecheck(self,col):
        rside = re.search("(?:[mM]erges?:)(.+)",col)
        if not rside:
            return {'withdrawn':False,'merges':[]}
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
            if (ind % 2) != 0: # before
                if rs == bfor:
                    correct = True
                    nxt.append(evenarr[ind])
        if correct:
            if aftr in nxt:
                after_match = True
        return {'correct':correct,'after_match':after_match,'alt_after':nxt}

    def mergedinto(self,evenarr,potbef):
        becomes = []
        for ind,rs in enumerate(evenarr):
            if (ind % 2) == 0: # before
                if potbef == rs:
                    becomes.append(evenarr[ind+1])
        return becomes

    def getb38(self,container):
        rslt = re.search('(?:b38=)(\S+)',container)
        if not rslt: # will need handling: AttributeError: 'NoneType' object has no attribute 'group'
            return None
        return rslt.group(1) 

    def add_alt(self,alt,main,ds):
        if self.reportmode:
            return
        vals = (main,alt)
        q = 'SELECT * FROM alt_ids WHERE id = %s AND alt_id = %s'
        self.omics.execute(q,vals)
        res = self.omics.fetchall()
        if len(res) > 0:
            res = ['|'.join(st) for st in res]
            self.dbact_om.write('-- failed to add %s as alt id to %s (ds %s) because it is already present: %s\n' % (alt,main,ds,'\t'.join(res)))
            return
        vals = (main,alt,ds)
        q = 'INSERT INTO alt_ids (id, alt_id, datasource) VALUES (%s,%s,%s)'
        try:
            self.omics.execute(q,vals)
        except:
            e = sys.exc_info()
            self.dbact_om.write( "-- Error trying to insert %s as alt id to %s\t%s\n" % (alt_id,main,e))

    def mtab_change_id(self,xsting,chngto,dbsnpin=False): # also for adding (vs replacing) new dbsnp id
        if xsting in self.actbr_idchng:
            print('entered mtab_change_id to change %s to %s but id %s has been edited already in this session' % (xsting,chngto,xsting))
            return
        if 'ukbbaffy_v2_1_map' in self.tabname:
            mid_col = 'chipid'
        else:
            mid_col = 'snp'
        where,val,colwithid,res = self.mtab_get_where_string(xsting)
        where = where.replace('%s','\'%s\'') + ';\n'
        if dbsnpin: # if unknown dbsnp id to be inserted instead of replaced (xsting mid maintained)
            colwithid = 'dbsnpid'
        if colwithid != mid_col:
            res = [rs[1] for rs in res]
            rs_ls = []
            for rs in res:
                rs_ls.extend(re.findall("rs[0-9]+", rs))
            rs_ls = list(dict.fromkeys(rs_ls))
            rs_ls = [rs for rs in rs_ls if rs != xsting and rs != chngto]
            rs_ls.append(chngto)
            chngto = ','.join(rs_ls)
        q = 'UPDATE ' + self.tabname + ' SET '  + colwithid + ' = \'' + chngto + '\'' + where
        self.dbact_br.write(q % val)
        self.actbr_idchng.add(xsting)

    def mtab_get_where_string(self,anid):
        if 'ukbbaffy_v2_1_map' in self.tabname:
            mid_col = 'chipid'
        else:
            mid_col = 'snp'
        rscol = 'dbsnpid'
        colwithid = ''
        where = ''
        res = None
        if 'rs' in anid: # or len(res)
            colwithid = rscol
            where = ' WHERE ' + rscol + ' REGEXP %s'
            q = 'SELECT ' + mid_col + ',dbsnpid,chr FROM ' + self.tabname + where
            val = (anid+'[[:>:]]',)
            self.br.execute(q,val)
            res = self.br.fetchall()
        if not res:
            colwithid = mid_col
            where = ' WHERE ' + colwithid + ' = %s'
            q = 'SELECT ' + mid_col + ',dbsnpid,chr FROM ' + self.tabname + where
            val = (anid,)
            self.br.execute(q,val)
            res = self.br.fetchall()
        if res:
            return where,val,colwithid,res
        else:
            return None

    def mtab_change_pos(self,anid,oldpos,newpos):
        if anid in self.actbr_poschng:
            print('entered mtab_change_pos but id %s has been edited already in this session' % (anid))
            return
        if oldpos in self.ignorepos or newpos in self.ignorepos:
            raise ValueError("ProcFile.mtab_change_pos(...) called to change map table %s position from %s to %s. Update to ignore %s needed" % (anid,oldpos,newpos,' and '.join(self.ignorepos)))
        where,val,rscol,res = self.mtab_get_where_string(anid)
        q = 'UPDATE ' + self.tabname + ' SET chr = %s ' + where + ' AND chr = %s'
        val = (newpos,) + val + (oldpos,)
        q = q.replace('%s','\'%s\'') + ';\n'
        self.dbact_br.write((q) % val)
        self.actbr_poschng.add(anid)
        #self.br.execute(q,val)
        #res = self.br.fetchall()
        #res = ['|'.join(row) for row in res]

    def swapout_main(self,swin,swout,ds):
        if self.reportmode:
            return True
        q = 'SELECT uid_datasource FROM consensus WHERE id = %s'
        val = (swout,)
        self.omics.execute(q,val)
        res = self.omics.fetchall()
        res = [rs[0] for rs in res if len(rs) > 0]
        if len(res) != 1:
            val = (swin,swout) + val
            self.dbact_om.write(('-- failed to swap in %s for %s because it may have been swapped already: ' + q + ';\n') % val)
            return False
        swout_ds = res[0]
        try:
            vals = (swin,ds,swout)
            q = 'UPDATE consensus SET id = %s, uid_datasource = %s WHERE id = %s'
            self.omics.execute(q,vals)
            vals = (swin,swout)
            q = 'UPDATE alt_ids SET id = %s WHERE id = %s'
            self.omics.execute(q,vals)
            q = 'UPDATE flank SET id = %s WHERE id = %s'
            self.omics.execute(q,vals)
            q = 'UPDATE positions SET id = %s WHERE id = %s'
            self.omics.execute(q,vals)
            q = 'UPDATE probes SET id = %s WHERE id = %s'
            self.omics.execute(q,vals)
            q = 'UPDATE snp_present SET id = %s WHERE id = %s'
            self.omics.execute(q,vals)
            q = 'UPDATE extra_map SET linkid = %s WHERE linkid = %s'
            self.omics.execute(q,vals)
            vals = (swin,swout,swout_ds)
            q = 'INSERT INTO alt_ids (id,alt_id,datasource) VALUES (%s,%s,%s)'
            self.omics.execute(q,vals)
        except:
            e = sys.exc_info()
            self.dbact_om.write( "Error trying to swap in %s for %s: %s:\n" % (swin,swout,e))
            self.swapout_main_write(swin,swout,ds,swout_ds)
            return False
        return True

    def swapout_main_write(self,swin,swout,ds,swout_ds):
        self.dbact_om.write('UPDATE consensus SET id = \'%s\', uid_datasource = \'%s\' WHERE id = \'%s\';\n' % (swin,ds,swout))
        self.dbact_om.write('UPDATE alt_ids SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE flank SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE positions SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE probes SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('UPDATE snp_present SET id = \'%s\' WHERE id = \'%s\';\n' % (swin,swout))
        self.dbact_om.write('INSERT INTO alt_ids (id, alt_id, datasource) VALUES (\'%s\',\'%s\',\'%s\');\n' % (swin,swout,swout_ds))
        self.dbact_om.write('--\n')

    def pos_flag(self,mid,chrpos,fl=-5):
        if self.reportmode:
            return
        vals = (fl,mid,'38',chrpos,-1)
        if chrpos in self.ignorepos:
            raise ValueError("FileProc.pos_flag(...) entered with chr:pos value %s" % (chrpos))
        q = 'UPDATE positions SET chosen = %s WHERE id = %s AND build = %s AND CONCAT(chr,":",pos) = %s AND chosen > %s'
        self.omics.execute(q,vals)

    def addpos(self,mid,chrpos,ds,build='38'):
        if self.reportmode:
            return
        if chrpos in self.ignorepos:
            raise ValueError("entered ProcFile.addpos(...) with chr:pos %s. Update required to ignore these strings" % (chrpos))
        q = 'SELECT * FROM positions WHERE id = %s AND build = %s AND CONCAT(chr,":",pos) = %s'
        vals = (mid,build,chrpos)
        self.omics.execute(q,vals)
        res = self.omics.fetchall()
        if len(res) > 0:
            res = ['|'.join([str(s) for s in st]) for st in res]
            self.dbact_om.write('-- Error adding position %s to %s (ds %s) because its already there: %s\n' % (chrpos,mid,ds,'\t'.join(res)))
            return
        chrm = chrpos.split(':')[0]
        pos = int(chrpos.split(':')[1])
        q = 'INSERT INTO positions (id, chr, pos, build, datasource) VALUES (%s,%s,%s,%s,%s)'
        vals = (mid,chrm,pos,build,ds)
        try:
            self.omics.execute(q,vals)
        except:
            e = sys.exc_info()
            self.dbact_om.write( "-- Error trying to add position %s to %s (ds %s): %s\n" % (chrpos,mid,ds,e))

    def checkbr_pos(self,mid,chrpos=None):
        allout = self.mtab_get_where_string(mid)
        if not allout:
            return None
        where,val,rscol,res = allout
        res = [s[2] for s in res] 
        res = list(dict.fromkeys(res))
        if not chrpos:
            return [None,res]
        if chrpos in res:
            return [True,res]
        return [False,res]

    def checkom_pos(self,mid,chrpos=None,build='38'):
        q = 'SELECT chr,pos FROM positions WHERE build = %s AND id = %s'
        val = (build,mid)
        self.omics.execute(q,val)
        res = self.omics.fetchall()
        res = [s[0] + ':'+ str(s[1]) for s in res]
        res = list(dict.fromkeys(res))
        if not chrpos:
            return [None,res]
        if chrpos in res:
            return [True,res]
        return [False,res]

    def checkomics_pos(self,mid):
        q = 'SELECT chr,pos FROM positions WHERE build = %s AND id = %s'
        val = ('38',mid)
        self.omics.execute(q,val)
        res = self.omics.fetchall()
        res = [s[0] + ':'+ str(s[1]) for s in res]
        return [res]

    def checkom_flag(self,mid,chrpos,build='38'):
        chrm,pos = chrpos.split(':')
        q = 'SELECT chosen FROM positions WHERE id = %s AND build = %s AND chr = %s AND pos = %s'
        vals = (mid,build,chrm,pos)
        self.omics.execute(q,vals)
        res = self.omics.fetchall()
        res = [ch[0] for ch in res]
        #res = [-1 for ch in res]
        badflag = False
        for chs in res:
            if chs < 0:
                badflag = True
        return [badflag,res]

    def checkomid(self,cid):
        q = 'SELECT EXISTS (SELECT 1 FROM consensus WHERE id = %s)'
        val = (cid,)
        self.omics.execute(q,val)
        conyesno = self.omics.fetchone()[0]
        q = 'SELECT EXISTS (SELECT 1 FROM alt_ids WHERE alt_id = %s)'
        self.omics.execute(q,val)
        altyesno = self.omics.fetchone()[0]
        return [conyesno,altyesno] # [known as main id, known at alternative id]

    def altomain(self,tobemain,oldmain):
        q = 'SELECT datasource FROM alt_ids WHERE alt_id = %s AND id = %s'
        vals = (tobemain,oldmain)
        self.omics.execute(q,vals)
        tobeds = [row[0] for row in self.omics.fetchall()]
        q = 'SELECT uid_datasource FROM consensus WHERE id = %s'
        vals = (oldmain,)
        self.omics.execute(q,vals)
        oldmds = [colo for colo in self.omics.fetchone()]
        print(tobeds,oldmds)
        if len(tobeds) != 1 or len(oldmds) != 1:
            return False
        #self.swapout_main(swin=tobemain,swout=oldmain,ds=tobeds)
        #delete from alt_ids where id = tobemain and alt_id = tobemain and datasource = tobeds
        return True

    def extra_map(self,newid,linkid,chrpos,datasource,chosen,ds_chrpos=None):
        if self.reportmode:
            return
        chrm = chrpos.split(':')[0]
        pos = chrpos.split(':')[1]
        q = "INSERT IGNORE INTO extra_map (newid,linkid,chr,GRCh38_pos,datasource,ds_chrpos,chosen) VALUES (%s,%s,%s,%s,%s,%s,%s)"
        vals = (newid,linkid,chrm,pos,datasource,ds_chrpos,chosen)
        if not ds_chrpos:
            q = "INSERT IGNORE INTO extra_map (newid,linkid,chr,GRCh38_pos,datasource,chosen) VALUES (%s,%s,%s,%s,%s,%s)"
            vals = (newid,linkid,chrm,pos,datasource,chosen)
        try:
            self.omics.execute(q,vals)
        except:
            e = sys.exc_info()
            self.dbact_om.write( "-- Error trying to add %s to extra_map under %s:\t%s\n" % (newid,linkid,e))

