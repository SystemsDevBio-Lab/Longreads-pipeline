#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import argparse

def main():
    ap = argparse.ArgumentParser(
        description="Create a gffutils SQLite DB from a (clean) GTF/GTF.gz"
    )
    ap.add_argument("gtf", help="Input GTF (can be .gtf or .gtf.gz)")
    ap.add_argument("--db", help="Output .db path (default: <gtf>.db)")
    ap.add_argument("--force", action="store_true", help="Overwrite existing DB")
    ap.add_argument("--merge", action="store_true",
                    help="Use merge_strategy=merge (by default use create_unique for safety)")
    ap.add_argument("--no-index", action="store_true",
                    help="Skip building extra helpful indexes (faster create, slower query)")
    args = ap.parse_args()
    
    gtf_file = args.gtf
    if not os.path.exists(gtf_file):
        print(f"[ERR] not found: {gtf_file}", file=sys.stderr)
        sys.exit(1)
    
    db_file = args.db or (gtf_file.rstrip(".gz") + ".db")
    
    if os.path.exists(db_file) and not args.force:
        print(f"[ERR] DB exists: {db_file}. Use --force to overwrite.", file=sys.stderr)
        sys.exit(2)
    
    try:
        import gffutils
    except ImportError:
        print("[ERR] Missing dependency: gffutils. Install with: pip install gffutils", file=sys.stderr)
        sys.exit(3)
    
    print(f"[INFO] creating DB from: {gtf_file}")
    print(f"[INFO] output DB: {db_file}")
    
    # 显式指定 ID 字段，与你清洗后的 GTF 对齐
    id_spec = {
        "gene": "gene_id",
        "transcript": "transcript_id",
        "exon": "exon_id",
    }
    
    merge_strategy = "merge" if args.merge else "create_unique"
    
    print(f"[INFO] merge_strategy: {merge_strategy}")
    
    db = gffutils.create_db(
        data=gtf_file,
        dbfn=db_file,
        force=args.force,
        keep_order=True,
        sort_attribute_values=True,
        merge_strategy=merge_strategy,
        id_spec=id_spec,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        # verbose=True,  # 需要更多日志时打开
    )
    
    if not args.no_index:
        print("[INFO] building helpful indexes ...")
        
        # 先检查数据库中有哪些表
        with db.conn as con:
            cursor = con.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tables = [row[0] for row in cursor.fetchall()]
            print(f"[DEBUG] Available tables: {tables}")
            
            # 基本索引（这些对所有版本都有效）
            try:
                con.execute("CREATE INDEX IF NOT EXISTS idx_features_seqid ON features(seqid)")
                print("[INFO] Created index: idx_features_seqid")
            except Exception as e:
                print(f"[WARN] Failed to create seqid index: {e}")
            
            try:
                con.execute("CREATE INDEX IF NOT EXISTS idx_features_type ON features(featuretype)")
                print("[INFO] Created index: idx_features_type")
            except Exception as e:
                print(f"[WARN] Failed to create featuretype index: {e}")
            
            try:
                con.execute("CREATE INDEX IF NOT EXISTS idx_features_start ON features(start)")
                print("[INFO] Created index: idx_features_start")
            except Exception as e:
                print(f"[WARN] Failed to create start index: {e}")
            
            try:
                con.execute("CREATE INDEX IF NOT EXISTS idx_features_end ON features(end)")
                print("[INFO] Created index: idx_features_end")
            except Exception as e:
                print(f"[WARN] Failed to create end index: {e}")
            
            # 尝试为 relations 表创建索引（如果存在）
            if 'relations' in tables:
                try:
                    con.execute("CREATE INDEX IF NOT EXISTS idx_relations_parent ON relations(parent)")
                    con.execute("CREATE INDEX IF NOT EXISTS idx_relations_child ON relations(child)")
                    print("[INFO] Created indexes for relations table")
                except Exception as e:
                    print(f"[WARN] Failed to create relations indexes: {e}")
            
            # gffutils 的 attributes 表名称可能是 'attributes' 或者属性存储在 features 表中
            # 检查 features 表的列
            cursor = con.execute("PRAGMA table_info(features)")
            columns = [row[1] for row in cursor.fetchall()]
            print(f"[DEBUG] Features table columns: {columns}")
            
            # 如果有 attributes 列（JSON格式），不需要额外的 attributes 表索引
            if 'attributes' in columns:
                print("[INFO] Attributes stored in features table (JSON format)")
            elif 'attributes' in tables:
                # 如果有单独的 attributes 表
                try:
                    con.execute("CREATE INDEX IF NOT EXISTS idx_attr_key ON attributes(key)")
                    con.execute("CREATE INDEX IF NOT EXISTS idx_attr_value ON attributes(value)")
                    print("[INFO] Created indexes for attributes table")
                except Exception as e:
                    print(f"[WARN] Failed to create attributes indexes: {e}")
    
    # 统计
    try:
        ng = db.count_features_of_type("gene")
    except:
        ng = 0
    
    try:
        nt = db.count_features_of_type("transcript")
    except:
        nt = 0
    
    try:
        ne = db.count_features_of_type("exon")
    except:
        ne = 0
    
    print(f"\n[DONE] {db_file}")
    print(f"[STAT] gene={ng}, transcript={nt}, exon={ne}")
    
    if not args.merge:
        print(f"[TIP] Default strategy is create_unique. If you truly need attribute merge for duplicated IDs, re-run with --merge.")

if __name__ == "__main__":
    main()
